#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/function.hpp"


using boost::lexical_cast;
namespace po=boost::program_options;

int main(int argc, char** argv){//main  

  //Input output and config options
  std::string cfg;
  bool concept;
  unsigned pNevts;
  std::string inFilePath;
  std::string outFilePath;
  unsigned nSiLayers;
  unsigned debug;

  po::options_description preconfig("Configuration"); 
  //config parameter input file
  //must contain one line per parameter with identifier=value pairs.
  preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    //Input output and config options //->required()
    ("concept",        po::value<bool>(&concept)->default_value(true))
    ("pNevts,n",       po::value<unsigned>(&pNevts)->default_value(0))
    ("inFilePath,i",   po::value<std::string>(&inFilePath)->required())
    ("outFilePath,o",  po::value<std::string>(&outFilePath)->required())
    ("nSiLayers",      po::value<unsigned>(&nSiLayers)->default_value(2))
    ("debug,d",        po::value<unsigned>(&debug)->default_value(0))
    ;

  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);

  //then run with 
  //./bin/executable -c configfile
  //or to run with a different option than the one in the config file:
  //./bin/executable -c configfile --inFilePath=/blah/blah/

  //each option has a short "-i", "-o" etc... 
  //specified in the lines above after the ",", 
  //and the long option with --inFilePath etc....



}
