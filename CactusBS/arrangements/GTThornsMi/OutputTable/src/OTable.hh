#ifndef OTABLE_IO_HEADER
#define OTABLE_IO_HEADER
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

class OutputTable
{
	private:
		int width;
		std::ofstream stream_obj;

	public:
		OutputTable(std::stringstream &fname, int w, std::ios_base::openmode md);

		~OutputTable();

		template<typename T>
			OutputTable& operator<<(const T& output);

		OutputTable& operator<<(std::ostream& (*func)(std::ostream&));
};


template<typename T>
OutputTable& OutputTable::operator<<(const T& output)
{
	stream_obj << std::setw(width) << output;
	return *this;	
}

#endif
