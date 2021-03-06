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

  double sumx=0.0;
  double sumy=0.0;
  double sumz=0.0;

  for(int i = 0; i < t; i++)
  {
    sumx += values.at(i).at(1);
    sumy += values.at(i).at(2);
    sumz += values.at(i).at(3);
  }

  double ave_x = sumx / t;
  double ave_y = sumy / t;
  double ave_z = sumz / t;

  //subtract off the average from each total
  //this reduces noise in the final output
  for(int i = 0; i < t; i++)
  {
    values.at(i).at(1) -= ave_x;
    values.at(i).at(2) -= ave_y;
    values.at(i).at(3) -= ave_z;
  }

  //window is the offset used in TCF.
  double window = 20000;

  //why an array? because doing the i-j business below
  //is a pain in the ass with vectors
  //(see the RHS of the assignments below)
  double * c = (double *)calloc(t, sizeof(double));

  for(int i = 0; i < t; i++)
  {
    if(i % 1000 == 0) { std::cout << "Step " << i << " reached." << std::endl; }
    for(int j = i; j < (i + window); j++)
    {
      if((i + window) > t) { break; }
      c[j-i] += (values.at(i).at(1) * values.at(j).at(1));
      c[j-i] += (values.at(i).at(2) * values.at(j).at(2));
      c[j-i] += (values.at(i).at(3) * values.at(j).at(3));
    }
  }

  std::string output_file_name = std::string(argv[1]) + ".CORRELATED";
  std::ofstream out(output_file_name);
  for(int i = 0; i < t; i++)
  {
    out << c[i]/(t-i) << std::endl;
  }

  free(c);

  return 0;
}
