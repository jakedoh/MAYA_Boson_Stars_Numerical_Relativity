#include <OutputTable.hh>

OutputTable::OutputTable(std::stringstream &fname, int w, std::ios_base::openmode md)
{
	stream_obj.open(fname.str().c_str(), md);
	width = w;
}

OutputTable::~OutputTable()
{
	stream_obj.close();
}

OutputTable& OutputTable::operator<<(std::ostream& (*func)(std::ostream&))
{
	func(stream_obj);
	return *this;
}
