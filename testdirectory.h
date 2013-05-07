#ifndef TESTDIRECTORY_H
#define TESTDIRECTORY_H
#include <iostream>
#include <dirent.h>

class TestDirectory
{
public:
    TestDirectory();
    bool DirectoryExists(const char *pzPath);
};

#endif // TESTDIRECTORY_H
