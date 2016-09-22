/***************************************************************************
 *  Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "kahypar/definitions.h"
#include "kahypar/io/hypergraph_io.h"
#include "kahypar/macros.h"

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cout << "No .hgr file specified" << std::endl;
  }
  std::string hgr_filename(argv[1]);
  std::string graphml_filename(hgr_filename + ".graph");

  Hypergraph hypergraph(partition::io::createHypergraphFromFile(hgr_filename, 2));

  std::ofstream out_stream(graphml_filename.c_str());

  // Vertices: hypernodes + hyperedges
  // Edges: One edge for each pin!

  out_stream << hypergraph.initialNumNodes() + hypergraph.initialNumEdges() << " " << hypergraph.initialNumPins() << std::endl;

  for (const HypernodeID hn : hypergraph.nodes()) {
    // vertex ids start with 1
    for (const HyperedgeID he : hypergraph.incidentEdges(hn)) {
      const HyperedgeID he_id = hypergraph.initialNumNodes() + he + 1;
      out_stream << he_id << " ";
    }
    out_stream << std::endl;
  }

  for (const HyperedgeID he : hypergraph.edges()) {
    for (const HypernodeID pin : hypergraph.pins(he)) {
      const HypernodeID hn_id = pin + 1;
      out_stream << hn_id << " ";
    }
    out_stream << std::endl;
  }

  out_stream.close();
  std::cout << " ... done!" << std::endl;
  return 0;
}
