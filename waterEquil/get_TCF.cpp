#include <iostream>
#include <fstream>
#include <vector>

int main(int argc, char ** argv)
{
  if(argc != 2)
  {
    std::cout << "Takes exactly one file name as input" << std::endl;
    exit(EXIT_FAILURE);
  }
  std::ifstream matrix(argv[1]);
  std::vector<std::vector<double>> values;
  std::vector<double> row(4);
  while(matrix >> row[0] >> row[1] >> row[2] >> row[3])
  {
    values.push_back(row);
  }

  int t = values.size();

  double * c = (double *)calloc(t, sizeof(double));

  for(int i = 0; i < t; i++)
  {
    for(int j = 0; j < i; j++)
    {
      c[i-j] += (values.at(i).at(1) * values.at(j).at(1));
      c[i-j] += (values.at(i).at(2) * values.at(j).at(2));
      c[i-j] += (values.at(i).at(3) * values.at(j).at(3));
    }
  }

  std::string output_file_name = std::string(argv[1]) + ".CORRELATED";
  std::ofstream out(output_file_name);
  for(int i = 0; i < t; i++)
  {
    out << " " <<  c[i]/(t-i) << std::endl;
  }

  return 0;
}
