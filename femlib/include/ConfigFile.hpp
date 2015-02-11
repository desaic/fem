#ifndef CONFIGFILE_HPP
#define CONFIGFILE_HPP
#include <map>
#include <string>
class ConfigFile
{
public:
  void load(const char * filename);
  float getFloat(const std::string & key);
  int getInt(const std::string & key);
  bool getBool(const std::string & key);
  std::string getString(const std::string & key);

  bool hasOpt(const std::string & key);

private:
  std::map<std::string, std::string> m;
};
#endif