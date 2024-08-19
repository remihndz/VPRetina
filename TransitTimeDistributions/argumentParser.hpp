#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <boost/lexical_cast.hpp>

void commandLineParser(int argc, char *argv[],		    
		       size_t &cutoff,
		       std::vector<std::string> &files){

  cutoff = 350; // Default value
  files.clear(); // Remove all elements
  
  std::string str;
  for (int i=1; i<argc; i++){
    str = argv[i];
    if ((str.find("-cutoff")!=std::string::npos)
	|| (str.find("-c")!=std::string::npos)){ // If (-)-cutoff or (-)-c in str
      i++; // Look at the argument following the flag
      cutoff = boost::lexical_cast<size_t>(argv[i]);
    }
    else
      files.push_back(str);
  }
};
