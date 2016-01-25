#include <windows.h>
#include<iostream>


#ifndef TIME_UTILITY_H
#define TIME_UTILITY_H

class TimeUtility{
private:

	double PCFreq;
	__int64 CounterStart;

public:
	TimeUtility();

	void StartCounterMicro();
	double GetCounterMicro();

	void StartCounterMill();
	double GetCounterMill();
};

#endif