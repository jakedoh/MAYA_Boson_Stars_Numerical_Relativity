#ifndef HTTPD_MAP_H_
#define HTTPD_MAP_H_ 1

#include <stddef.h>

#if __cplusplus
extern "C" {
#define VOID
#else
#define VOID void
#endif

typedef void * uMap;

uMap Httpd_MapCreate(VOID);
void * Httpd_MapData(uMap map, size_t keylen, const char * key);
int Httpd_MapStore(uMap map, size_t keylen, const char * key, void * data);
void Httpd_MapDestroy(uMap map, void (*destroy)(void *));

#if __cplusplus
}
#endif

#endif
