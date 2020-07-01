/*
 * OctoMap - An Efficient Probabilistic 3D Mapping Framework Based on Octrees
 * http://octomap.github.com/
 *
 * Copyright (c) 2009-2013, K.M. Wurm and A. Hornung, University of Freiburg
 * All rights reserved.
 * License: New BSD
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Freiburg nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <octomap/VirtualOcTree.h>

namespace octomap {

// node implementation  --------------------------------------
std::ostream& VirtualOcTreeNode::writeData(std::ostream& s) const {
  s.write((const char*)&value, sizeof(value));      // occupancy
  s.write((const char*)&is_real, sizeof(is_real));  // is_real

  return s;
}

std::istream& VirtualOcTreeNode::readData(std::istream& s) {
  s.read((char*)&value, sizeof(value));      // occupancy
  s.read((char*)&is_real, sizeof(is_real));  // is_real

  return s;
}

bool VirtualOcTreeNode::getMajorityChildIsReal() const {
  uint8_t count_is_real = 0u;

  if (children != NULL) {
    for (unsigned int i = 0; i < 8; i++) {
      if (children[i] != NULL) {
        bool is_real =
            static_cast<VirtualOcTreeNode*>(children[i])->getIsReal();
        if (is_real) {
          ++count_is_real;
        }
      }
    }
  }

  return count_is_real > 4u;
}

VirtualOcTree::VirtualOcTree(double in_resolution)
    : OccupancyOcTreeBase<VirtualOcTreeNode>(in_resolution) {
  virtualOcTreeMemberInit.ensureLinking();
};

VirtualOcTree::StaticMemberInitializer VirtualOcTree::virtualOcTreeMemberInit;

void VirtualOcTree::insertPointCloud(const Pointcloud& scan,
                                     const octomap::point3d& sensor_origin,
                                     bool is_real, double maxrange,
                                     bool lazy_eval, bool discretize) {
  KeySet free_cells, occupied_cells;
  if (discretize)
    computeDiscreteUpdate(scan, sensor_origin, free_cells, occupied_cells,
                          maxrange);
  else
    computeUpdate(scan, sensor_origin, free_cells, occupied_cells, maxrange);
  bool occupied = false;

  // insert data into tree  -----------------------
  for (KeySet::iterator it = free_cells.begin(); it != free_cells.end(); ++it) {
    updateNode(*it, occupied, is_real, lazy_eval);
  }

  occupied = true;
  for (KeySet::iterator it = occupied_cells.begin(); it != occupied_cells.end();
       ++it) {
    updateNode(*it, occupied, is_real, lazy_eval);
  }
}

VirtualOcTreeNode* VirtualOcTree::updateNode(const OcTreeKey& key,
                                             bool occupied, bool is_real,
                                             bool lazy_eval) {
  float logOdds = this->prob_miss_log;
  if (occupied) logOdds = this->prob_hit_log;

  return updateNode(key, logOdds, is_real, lazy_eval);
}

VirtualOcTreeNode* VirtualOcTree::updateNode(const OcTreeKey& key,
                                             float log_odds_update,
                                             bool is_real, bool lazy_eval) {
  // TODO(margaritaG): disabled.
  // early abort (no change will happen).
  // may cause an overhead in some configuration, but more often helps
  // VirtualOcTreeNode* leaf = this->search(key);
  // // no change: node already at threshold
  // if (leaf && ((log_odds_update >= 0 &&
  //               leaf->getLogOdds() >= this->clamping_thres_max) ||
  //              (log_odds_update <= 0 &&
  //               leaf->getLogOdds() <= this->clamping_thres_min))) {
  //   return leaf;
  // }

  bool createdRoot = false;
  if (this->root == NULL) {
    this->root = new VirtualOcTreeNode();
    this->tree_size++;
    createdRoot = true;
  }

  return updateNodeRecurs(this->root, createdRoot, key, 0, log_odds_update,
                          is_real, lazy_eval);
}

VirtualOcTreeNode* VirtualOcTree::updateNodeRecurs(
    VirtualOcTreeNode* node, bool node_just_created, const OcTreeKey& key,
    unsigned int depth, const float& log_odds_update, bool is_real,
    bool lazy_eval) {
  bool created_node = false;

  assert(node);

  // follow down to last level
  if (depth < this->tree_depth) {
    unsigned int pos = computeChildIdx(key, this->tree_depth - 1 - depth);
    if (!this->nodeChildExists(node, pos)) {
      // child does not exist, but maybe it's a pruned node?
      if (!this->nodeHasChildren(node) && !node_just_created) {
        // current node does not have children AND it is not a new node
        // -> expand pruned node
        this->expandNode(node);
      } else {
        // not a pruned node, create requested child
        this->createNodeChild(node, pos);
        created_node = true;
      }
    }

    if (lazy_eval)
      return updateNodeRecurs(this->getNodeChild(node, pos), created_node, key,
                              depth + 1, log_odds_update, is_real, lazy_eval);
    else {
      VirtualOcTreeNode* retval =
          updateNodeRecurs(this->getNodeChild(node, pos), created_node, key,
                           depth + 1, log_odds_update, is_real, lazy_eval);
      // prune node if possible, otherwise set own probability
      // note: combining both did not lead to a speedup!
      if (this->pruneNode(node)) {
        // return pointer to current parent (pruned), the just updated node no
        // longer exists
        retval = node;
      } else {
        node->updateOccupancyChildren();
      }

      return retval;
    }
  }

  // at last level, update node, end of recursion
  else {
    if (use_change_detection) {
      bool occBefore = this->isNodeOccupied(node);
      updateNodeLogOdds(node, log_odds_update);

      if (node_just_created) {  // new node
        changed_keys.insert(std::pair<OcTreeKey, bool>(key, true));
      } else if (occBefore !=
                 this->isNodeOccupied(node)) {  // occupancy changed, track it
        KeyBoolMap::iterator it = changed_keys.find(key);
        if (it == changed_keys.end())
          changed_keys.insert(std::pair<OcTreeKey, bool>(key, false));
        else if (it->second == false)
          changed_keys.erase(it);
      }
    } else {
      updateNodeLogOdds(node, log_odds_update);
      updateNodeIsReal(node, is_real);
    }
    return node;
  }
}

void VirtualOcTree::updateNodeIsReal(VirtualOcTreeNode* virtualOccupancyNode,
                                     const bool& update) const {
  // TODO(margaritaG): is this the best way?
  // For now simple strategy: only update the node is_real value to true if
  // previously it was false and the new update (is_real) has value true.
  if (!virtualOccupancyNode->getIsReal() && update) {
    virtualOccupancyNode->setIsReal(update);
  }
}

}  // namespace octomap
