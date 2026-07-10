#ifndef LISTHELPERS_H
#define LISTHELPERS_H

#include "../common/defines.h"
#include <vector>

template<typename T>
static inline bool AnyNullPtr(std::vector<T> list)
{
    for(int i = 0; i < list.size(); i++)
    {
        if(list[i] == nullptr)
        {
            return true;
        }
    }
    return false;
}
template<typename T>
static inline void SafeDeleteList(std::vector<T> list)
{
    for(int i = 0; i < list.size(); i++)
    {
        SAFE_DELETE(list[i]);
    }
}

#endif // LISTHELPERS_H
