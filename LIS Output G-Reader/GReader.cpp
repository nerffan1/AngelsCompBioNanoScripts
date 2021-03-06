#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

//Description: This program reads the G Section of a .lis file generated by Curves5.3
//Argument: Name of file to read
//Example: g++ -c GReader.cpp -o THIS_IS_A_FILE.lis (?)
//LAST UPDATED: 03/22/20 at 9:47pm

//Constants
const size_t STRAND_NUM = 4;
const size_t PAIR_NUM = 2;
const size_t PARAM_NUM = 6;

//Additional Functions
void PrintToCSV(std::stringstream &Line, std::fstream & MyCSV);

int main(int argc, char* argv[])
{
	//Open Files
	std::fstream MyLisFile, MyCSVFile;
	MyLisFile.open(argv[1], std::ios::in);
	MyCSVFile.open(argv[2], std::ios::out | std::ios::app);

	//Read information
	if (MyLisFile.is_open())
	{
		//Find the section with |G|
		std::string Line;
		while (Line.find("|G|") == std::string::npos)
		{
			std::getline(MyLisFile, Line);
		}
		std::getline(MyLisFile, Line); //Jump one more line to start at the "------" line

		//Loop through the 4 strands
		for (size_t i = 0; i < STRAND_NUM; i++)
		{
			//Jump to line before info
			for (size_t i = 0; i < 4; i++)
			{
				std::getline(MyLisFile, Line);
			}
			//Read info lines
			for (size_t i = 0; i < PAIR_NUM; i++)
			{
				std::getline(MyLisFile, Line);
				std::stringstream MyString;
				MyString << Line.substr(Line.find_first_of("   ",16)); //This number was chosen arbitrarily
				std::cout << MyString.str() << std::endl;
				PrintToCSV(MyString,MyCSVFile);
			}
		}
	}
	MyLisFile.close();
	MyCSVFile.close();
	return 0;
}

void PrintToCSV(std::stringstream &Line, std::fstream &MyCSV)
{
	for (size_t i = 0; i < PARAM_NUM; i++)
	{
		float num;
		Line >> num;
		MyCSV << num << ",";
	}
	MyCSV << std::endl;
}
