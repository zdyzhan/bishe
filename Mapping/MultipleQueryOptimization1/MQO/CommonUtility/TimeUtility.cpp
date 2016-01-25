#include"TimeUtility.h"




TimeUtility::TimeUtility()
{
	PCFreq = 0.0;
	CounterStart = 0;
}

void TimeUtility::StartCounterMicro()
{
    LARGE_INTEGER li;
    if(QueryPerformanceFrequency(&li))
		PCFreq = li.QuadPart/1000000.0;
    QueryPerformanceCounter(&li);
    CounterStart = li.QuadPart;
}
double TimeUtility::GetCounterMicro()
{
    LARGE_INTEGER li;
    QueryPerformanceCounter(&li);
    return double(li.QuadPart-CounterStart)/PCFreq;
}

void TimeUtility::StartCounterMill(){
    LARGE_INTEGER li;
    if(QueryPerformanceFrequency(&li))
		PCFreq = li.QuadPart/1000.0;
    QueryPerformanceCounter(&li);
    CounterStart = li.QuadPart;
}

double TimeUtility::GetCounterMill(){
    LARGE_INTEGER li;
    QueryPerformanceCounter(&li);
    return double(li.QuadPart-CounterStart)/PCFreq;
}