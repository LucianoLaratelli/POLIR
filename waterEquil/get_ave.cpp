#include <iostream>
#include <vector>
#include <fstream>


int main(int argc, char ** argv)
{
  if (argc != 2) { exit(EXIT_FAILURE); }
  std::ifstream in_file(argv[1]);
  int step;
  double energy;
  double energy_sum = 0.0;
  std::vector <double> average;

  std::ofstream out_file("energy_ave.txt");

  while(in_file >> step >> energy)
  {
    energy_sum += energy;
    if(step == 0) { continue; }
    else
    {
      out_file << step  << " "<< energy_sum / step << std::endl;
    }
  }

  in_file.close();
  out_file.close();

  return 0;
}
