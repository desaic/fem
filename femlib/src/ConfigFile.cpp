#include "ConfigFile.hpp"
#include "FileUtil.hpp"
#include <sstream>

void ConfigFile::load(const char * filename)
{
  FileUtilIn in(filename);
  if (!in.good()){
    return;
  }
  std::string line;
  while (std::getline(in.in, line)){
    std::istringstream ss(line);
    if (line.size() < 3){
      continue;
    }
    std::string token;
    ss >> token;
    if (token[0] == '#'){
      continue;
    }
    std::string val;
    ss >> val;
    m[token] = val;
    std::cout << token << ": " << m[token] << "\n";
  }
  in.close();
}

float ConfigFile::getFloat(const std::string & key)
{
  if (!hasOpt(key)){
    return 0.0f;
  }
  return std::stof(m[key]);
}

int ConfigFile::getInt(const std::string & key)
{
  if (!hasOpt(key)){
    return 0;
  }
  return std::stoi(m[key]);
}

std::string ConfigFile::getString(const std::string & key)
{
  if (!hasOpt(key)){
    return std::string();
  }
  return m[key];
}

bool ConfigFile::getBool(const std::string & key)
{
  if (!hasOpt(key)){
    return false;
  }
  return m[key] == "true";
}

bool ConfigFile::hasOpt(const std::string & key)
{
  return m.find(key) != m.end();
}
