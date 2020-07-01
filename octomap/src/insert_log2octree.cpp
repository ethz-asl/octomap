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
#include <octomap/octomap.h>
#include <stdlib.h>
#include <string.h>

using namespace std;
using namespace octomap;

void printUsage(char* self) {
  std::cerr << "\nUSAGE: " << self << " [options]\n\n";

  std::cerr << "This tool inserts the point cloud described by a plain text "
               "log file into a virtual octree."
            << std::endl;
  std::cerr << "The log file needs to be in the format of:\n"
            << "NODE x y z roll pitch yaw\n"
            << "x y z\nx y z\n...\n"
            << "NODE x y z roll pitch yaw\n"
            << "x y z\n...\n\n"
            << "Lines starting with '#' or empty lines are ignored.\n\n";

  std::cerr << "OPTIONS:\n  -i <InputTree.bt> (required)\n"
               "  -log <InputCloud.log> (required) \n"
               "  -o <OutputTree.bt> (required) \n"
               "  -m <maxrange> (optional) \n"
               "  -is_predicted (optional, default: false) \n"
               "\n";

  exit(0);
}

int main(int argc, char** argv) {
  // default values:
  string inputTreeFilename = "";
  string logFilename = "";
  string outputTreeFilename = "";
  double maxrange = -1;
  bool discretize = false;
  bool lazy_eval = false;
  bool is_real = true;

  int arg = 0;
  while (++arg < argc) {
    if (!strcmp(argv[arg], "-i"))
      inputTreeFilename = std::string(argv[++arg]);
    else if (!strcmp(argv[arg], "-log"))
      logFilename = std::string(argv[++arg]);
    else if (!strcmp(argv[arg], "-o"))
      outputTreeFilename = std::string(argv[++arg]);
    else if (!strcmp(argv[arg], "m") && argc - arg < 3)
      printUsage(argv[0]);
    else if (!strcmp(argv[arg], "-m"))
      maxrange = atof(argv[++arg]);
    else if (!strcmp(argv[arg], "-is_predicted"))
      is_real = false;
    else {
      printUsage(argv[0]);
    }
  }

  cout << "\nReading Tree file\n===========================\n";

  AbstractOcTree* readTreeAbstract = AbstractOcTree::read(inputTreeFilename);
  VirtualOcTree* tree = dynamic_cast<VirtualOcTree*>(readTreeAbstract);

  cout << "\nReading Log file\n===========================\n";
  ScanGraph* graph = new ScanGraph();
  graph->readPlainASCII(logFilename);

  cout << "\nInserting point cloud into "
          "tree\n===========================\n";

  size_t numScans = graph->size();
  unsigned int currentScan = 1;
  for (ScanGraph::iterator scan_it = graph->begin(); scan_it != graph->end();
       scan_it++) {
    tree->insertPointCloud((*scan_it)->scan, (*scan_it)->pose.trans(), is_real,
                          maxrange, lazy_eval, discretize);
    cout << "(" << currentScan << "/" << numScans << ") " << flush;

    currentScan++;
  }

  // get rid of graph in mem before doing anything fancy with tree (=> memory)
  delete graph;

  cout << "\nDone inserting into tree.\n\n";

  cout << "\nWriting tree file\n===========================\n";
  tree->write(outputTreeFilename);

  return 0;
}
