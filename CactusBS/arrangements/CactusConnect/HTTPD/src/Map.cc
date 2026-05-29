#include <map>
#include <string>

#include "httpd_Map.h"

typedef std::map<std::string,void *> iMap;

extern "C" {

uMap Httpd_MapCreate(void)
{
  return new iMap;
}

void * Httpd_MapData(uMap map, size_t keylen, const char * key)
{
  iMap& imap(*static_cast<iMap*>(map));
  const std::string ikey(key, keylen);
  iMap::iterator it(imap.find(ikey));

  if(it != imap.end())
    return it->second;
  else
    return NULL;
}

int Httpd_MapStore(uMap map, size_t keylen, const char * key, void * data)
{
  iMap& imap(*static_cast<iMap*>(map));
  const std::string ikey(key, keylen);
  iMap::iterator it(imap.find(ikey));

  if(it != imap.end())
    it->second = data;
  else
    imap[ikey] = data;

  return 0;
}

void Httpd_MapDestroy(uMap map, void (*destroy)(void *))
{
  iMap& imap(*static_cast<iMap*>(map));
  for(iMap::iterator it = imap.begin() ; it != imap.end() ; ++it)
    destroy(it->second);
  delete &imap;
}

}
