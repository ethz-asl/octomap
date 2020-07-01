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

#ifndef OCTOMAP_VIRTUAL_OCTREE_H
#define OCTOMAP_VIRTUAL_OCTREE_H

#include <octomap/OcTreeNode.h>
#include <octomap/OccupancyOcTreeBase.h>
#include <iostream>

namespace octomap {

// forward declaraton for "friend"
class VirtualOcTree;

// node definition
class VirtualOcTreeNode : public OcTreeNode {
 public:
  friend class VirtualOcTree;  // needs access to node children (inherited)

 public:
  VirtualOcTreeNode() : OcTreeNode(), is_real(false) {}

  bool operator==(const VirtualOcTreeNode& rhs) const {
    return (rhs.value == value && rhs.is_real == is_real);
  }

  void copyData(const VirtualOcTreeNode& from) {
    OcTreeNode::copyData(from);
    this->is_real = from.getIsReal();
  }

  inline void setIsReal(bool is_r) { is_real = is_r; }

  inline bool getIsReal() const { return is_real; }

  /**
   * @return maximum of children's occupancy probabilities, in log odds
   */
  bool getMajorityChildIsReal() const;

  /// update this node's occupancy according to its children's maximum occupancy
  inline void updateOccupancyChildren() {
    this->setLogOdds(this->getMaxChildLogOdds());  // conservative
    // TODO(margaritaG): is this the right thing to do?
    // Maybe to be conservative we should set is_real to false if even one of
    // the children is not real?
    this->setIsReal(this->getMajorityChildIsReal());
  }

  // file I/O
  std::istream& readData(std::istream& s);
  std::ostream& writeData(std::ostream& s) const;

 protected:
  bool is_real;
};

// tree definition
class VirtualOcTree : public OccupancyOcTreeBase<VirtualOcTreeNode> {
 public:
  /// Default constructor, sets resolution of leafs
  VirtualOcTree(double resolution);

  /// virtual constructor: creates a new object of same type
  /// (Covariant return type requires an up-to-date compiler)
  VirtualOcTree* create() const { return new VirtualOcTree(resolution); }

  std::string getTreeType() const { return "VirtualOcTree"; }

  /**
   * Integrate a Pointcloud (in global reference frame), parallelized with
   * OpenMP. Special care is taken that each voxel in the map is updated only
   * once, and occupied nodes have a preference over free ones. This avoids
   * holes in the floor from mutual deletion and is more efficient than the
   * plain ray insertion in insertPointCloudRays().
   *
   * @note replaces insertScan()
   *
   * @param scan Pointcloud (measurement endpoints), in global reference frame
   */
  virtual void insertPointCloud(const Pointcloud& scan,
                                const octomap::point3d& sensor_origin,
                                bool is_real, double maxrange, bool lazy_eval,
                                bool discretize);

  /**
   * Integrate occupancy measurement.
   *
   * @param key OcTreeKey of the NODE that is to be updated
   * @param occupied true if the node was measured occupied, else false
   * @param lazy_eval whether update of inner nodes is omitted after the update
   * (default: false). This speeds up the insertion, but you need to call
   * updateInnerOccupancy() when done.
   * @return pointer to the updated NODE
   */
  virtual VirtualOcTreeNode* updateNode(const OcTreeKey& key, bool occupied,
                                        bool is_real, bool lazy_eval = false);

  virtual VirtualOcTreeNode* updateNode(const OcTreeKey& key,
                                        float log_odds_update, bool is_real,
                                        bool lazy_eval = false);

  /// update is_real value of node.
  virtual void updateNodeIsReal(VirtualOcTreeNode* virtualOccupancyNode,
                                const bool& update) const;

 protected:
  // recursive calls ----------------------------

  VirtualOcTreeNode* updateNodeRecurs(VirtualOcTreeNode* node,
                                      bool node_just_created,
                                      const OcTreeKey& key, unsigned int depth,
                                      const float& log_odds_update,
                                      bool is_real, bool lazy_eval = false);

  /**
   * Static member object which ensures that this OcTree's prototype
   * ends up in the classIDMapping only once. You need this as a
   * static member in any derived octree class in order to read .ot
   * files through the AbstractOcTree factory. You should also call
   * ensureLinking() once from the constructor.
   */
  class StaticMemberInitializer {
   public:
    StaticMemberInitializer() {
      VirtualOcTree* tree = new VirtualOcTree(0.1);
      tree->clearKeyRays();
      AbstractOcTree::registerTreeType(tree);
    }

    /**
     * Dummy function to ensure that MSVC does not drop the
     * StaticMemberInitializer, causing this tree failing to register.
     * Needs to be called from the constructor of this octree.
     */
    void ensureLinking(){};
  };

  /// static member to ensure static initialization (only once)
  static StaticMemberInitializer virtualOcTreeMemberInit;
};

}  // namespace octomap

#endif
