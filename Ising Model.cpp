/*
This program uses the Metropolis-Hastings algorithm to find numerical data for the 2D ising model on a square lattice.  Quantities that can be found with this program are: Average magnetization, average energy density, magnetic suscpetibility, specific heat at different values of J,h and T and for various lattice sizes.

The code also can be used to demonstrate the phenomenon of hysteresis.

Samples of around 100000 should be used to get good data

Equilibration requires between 500 and 1000 initial sweeps of the Metropolis-Hastings algorithm
*/

#include <stdio.h>
#include<stdlib.h>
#include<random>
#include<vector>
#include<cmath>
#include<iostream>
#include <time.h>
#include <fstream>
#include <numeric>

std::random_device rd;
std::mt19937 gen(rd());

std::vector<std::vector<int>> generate_lattice(int lattice_size)
{
    //Generates a random initial 2D lattice for monte carlo
    std::vector<std::vector<int>> lattice(lattice_size, std::vector<int>(lattice_size));
    std::uniform_int_distribution<> spin_choice(0,1);//Used to randomly assign spin values to each lattice site
    for (int i = 0; i < lattice.size(); ++i)//Initialise lattice in some random configuration
    {
        for (int j = 0; j< lattice[i].size(); ++j)
        {
            lattice[i][j] = pow(-1,spin_choice(gen));
        }
    }
    return lattice;
}

std::vector<float> linspace(float a, float b, int num)
{
             // create a vector of length num
             std::vector<float> v(num);
             double tmp = 0.0;
             
             // now assign the values to the vector
             for (int i = 0; i < num; i++)
                 {
                      v[i] = a + i * (float(b - a) / num);
                 }
             return v;
}
    
class Lattice
{
    
    public:
    int lattice_size; //Length of lattice, number of sites = lattice_size^2
    float T; //Temperature variable
    float J; //Site coupling variable
    float h; //Magnetic field variable
    int warm_up_iterations; //Number of iterations to bring random initial state to equilibrium
    int samples; //Number of monte carlo samples for which data is collected
    std::vector<std::vector<int>> lattice;
    
    Lattice(int size_of_lattice, float temperature, float coupling, float magnetic_field, int equilibration_iterations, int mc_samples, std::vector<std::vector<int>> initial_lattice)//Class initialization constructor
    {
        lattice_size = size_of_lattice;
        T = temperature; 
        J = coupling; 
        h = magnetic_field;
        warm_up_iterations = equilibration_iterations;
        samples =mc_samples;
        lattice = initial_lattice;
    }
    
    int neighbours(int x_position,int y_position)
    {
        //Finds the total spin of nearest neighbours of site (x_position,y_position)
        int R = lattice[(x_position+1)%lattice.size()][y_position];
        int L = lattice[(x_position-1)%lattice.size()][y_position];
        int U = lattice[x_position][(y_position+1)%lattice.size()];
        int D = lattice[x_position][(y_position-1)%lattice.size()];
        
        return R+L+U+D;
    }

    float dE(int x_position,int y_position)
    {
        //Finds change in energy from flipping the spin at site (x_position,y_position)
        return -2*lattice[x_position][y_position]*neighbours(x_position,y_position)*J - 2*h*lattice[x_position][y_position];
    }

    void flips()
    {
        //Goes through the lattice and randomly flips N^2 spins for a NxN lattice
        std::uniform_real_distribution<double> flip_choice(0,1);
        std::uniform_int_distribution<> site_choice(0,lattice.size()-1);
        
        for(int flip = 0; flip < pow(lattice.size(),2);++flip)
        {
            int a = site_choice(gen);//Choose random site and flip it
            int b = site_choice(gen);
            lattice[a][b] = -lattice[a][b];
            float delta_E = dE(a,b);//find energy and accept if dE < 0 s.t the flip is favourable
            if(delta_E > 0)//otherwise only accept the flip with some probability weighted by the temperature
            {
                float p = std::exp(-delta_E/(T));
                if(flip_choice(gen)>p)
                {
                    lattice[a][b] = -lattice[a][b];
                }
            }
        }
    }
    
    std::vector<std::vector<float>> get_energy_and_magnetization()
    {
        //Finds total absolute magnetization and energy for a given spin configuration for each monte carlo sample
        //Found simultaneously to avoid redundantly running monte Carlo twice
        std::vector<float> magnetization(samples);
        std::vector<float> energy(samples);
        std::vector<std::vector<float>> energy_and_magnetization(2, std::vector<float>(samples));
        for(int sample = 0; sample<samples; ++sample)
        {
            float total_spin = 0;
            float coupling_term = 0;
            float onsite_term = 0;
            for(int x_position = 0; x_position<lattice.size();++x_position)
            {
                for(int y_position = 0; y_position<lattice.size();++y_position)
                {
                    coupling_term += lattice[x_position][y_position]*neighbours(x_position,y_position);
                    onsite_term += lattice[x_position][y_position];
                    total_spin += lattice[x_position][y_position];
                }
            }
            magnetization[sample] = abs(total_spin);
            energy[sample] = -J*coupling_term/2 + h*onsite_term;//factor of 1/2 takes into account double counting of sites for coupling term
            flips();//after each iteration update the state by performing metropolis algorithm on L^2 random spins
        }
        energy_and_magnetization[0] = energy;
        energy_and_magnetization[1] = magnetization;
        
        return energy_and_magnetization;
    }
    
    float get_susceptibility(std::vector<float> magnetization)
    {   
        //Finds magnetic suscpetibility which can be shown to be Var(M)/T.  
        //Variance is with respect to different spin configurations.
        //magnetization holds total magnetization for some number of different equilibrated spin configurations
        //the following finds variance of that array
        float mag_avg = 0;
        float variance = 0;
        
        for(int i = 0; i < magnetization.size(); ++ i)
        {
            mag_avg+=magnetization[i]/magnetization.size();
        }
        
        for(int i = 0; i<magnetization.size(); ++i)
        {
            variance += pow((magnetization[i]-mag_avg),2)/magnetization.size();
        }
        
        return variance/T;
    }

    float get_specific_heat(std::vector<float> energy)//Not working for some reason
    {
        //Calculates specific heat at equilibrium which is given by var(E)/T^2
        //Variance is with respect to different spin configurations.
        //energy holds total energy for some number of different equilibrated spin configurations
        //the following finds variance of that array
        float avg_energy = 0;
        float variance = 0;
        
        for(int i = 0; i < energy.size(); ++i)
        {
            avg_energy += energy[i]/energy.size();
        }
        
        for(int i = 0; i < energy.size(); ++i)
        {
            variance += pow(energy[i]-avg_energy,2)/energy.size();
        }
        
        return variance/pow(T,2);
        
    }
    
    void warm_up_lattice()
    {
        //brings initial random lattice to equilibrium by performing a large number of metropolis sweeps
        for(int sweep_number = 0; sweep_number<warm_up_iterations; ++sweep_number)
        {
            flips();
        }
    }
    
    void print_lattice()
    {
        //Just prints current lattice configuration.  Used for testing and debugging
        for(int i = 0; i<lattice_size; ++i)
        {
            for(int j = 0; j<lattice_size;++j)
            {
                std::cout<<lattice[i][j]<<" ";
            }
            std::cout<<"\n";
        }
    }
    
    std::vector<float> get_hysteresis_loop(std::vector<float> magnetic_field)
    {
        std::vector<float> hysteresis_loop(magnetic_field.size());
        for(int h_index = 0; h_index < magnetic_field.size(); ++h_index)
        {
            h = magnetic_field[h_index];
            flips();
            for(int x_position = 0; x_position<lattice.size();++x_position)
            {
                for(int y_position = 0; y_position<lattice.size();++y_position)
                {
                    hysteresis_loop[h_index] += lattice[x_position][y_position]/pow(lattice.size(),2);
                }
            }   
        }
        return hysteresis_loop;
    }
};

int main()
{
    clock_t tStart = clock();
    
    int size_data_points = 8; //Number of lattice sizes that are looked at
    int temp_data_points = 100; // number of data points for whatever quantities are being found
    
    float T_init = 1.5;//Starting temperature
    float T_fin  = 3; //Final temperature
    float dT = (T_fin - T_init)/temp_data_points;//Temperature increment

    //Below stores quantities for different lattice sizes and temperatures
    std::vector<std::vector<float>> susceptibility(size_data_points, std::vector<float>(temp_data_points));
    std::vector<std::vector<float>> magnetization(size_data_points, std::vector<float>(temp_data_points));
    std::vector<std::vector<float>> energy(size_data_points, std::vector<float>(temp_data_points));
    std::vector<std::vector<float>> specific_heat(size_data_points, std::vector<float>(temp_data_points));
    
    int lattice_size = 10;//Number of sites is lattice_size^2
    std::vector<std::vector<int>> initial_lattice(lattice_size, std::vector<int>(lattice_size));

    float d = 0;
    for(int size_index = 0;size_index<size_data_points;++size_index)
    { 
        float T = T_init;
        //following loop finds per site magnetization, energy density, specific heat and magnetic suscpetibility for various temperaturs
        for(int temp_index = 0;temp_index<temp_data_points;++temp_index)
        {
            initial_lattice = generate_lattice(lattice_size);//Initialise lattice
            Lattice obj1(lattice_size,T , 1, 0, 750,100000,initial_lattice); //Initialise Ising model class
            obj1.warm_up_lattice();//Equilibrate lattice with "warm_up_iterations" number of sweeps of M-H algorithm
            
            std::vector<std::vector<float>> temp = obj1.get_energy_and_magnetization();//Stores energy and magnetization
            std::vector<float> temporary_energy = temp[0];//Stores the energy for each MC sample for this specific temperature and size
            std::vector<float> temporary_magnetization = temp[1];//Stores the magnetization for each MC sample for this specific temperature and size
            
            
            susceptibility[size_index][temp_index] = obj1.get_susceptibility(temporary_magnetization)/pow(lattice_size,2);
            magnetization[size_index][temp_index]  = std::accumulate(temporary_magnetization.begin(), temporary_magnetization.end(),0.0)/(temporary_magnetization.size()*pow(lattice_size,2));
            
            specific_heat[size_index][temp_index] = obj1.get_specific_heat(temporary_energy)/pow(lattice_size,2);
            energy[size_index][temp_index]  = std::accumulate(temporary_energy.begin(), temporary_energy.end(),0.0)/(temporary_energy.size()*pow(lattice_size,2));
            
            T+=dT;
            std::cout<<float(d)*100/(temp_data_points*size_data_points)<<"%\r";
            std::cout.flush();
            d++;
        }
        lattice_size+=2;
    }
//Uncomment below to generate a hysteresis loop
/*    
    the following generates a hysteresis loop at temperature T between the range(-h_fin,h_fin)
    int data_points = 30; // 1/3 of number of data points for hysteresis loop 

    T = 2;
    float h_fin = 6.5;
    std::vector<float> magnetic_field(linspace(0,h_fin,data_points));
    std::vector<float> magnetic_field2(linspace(h_fin,-h_fin,data_points));
    std::vector<float> magnetic_field3(linspace(-h_fin,h_fin,data_points));
    
    magnetic_field.insert(magnetic_field.end(),magnetic_field2.begin(),magnetic_field2.end());
    magnetic_field.insert(magnetic_field.end(),magnetic_field3.begin(),magnetic_field3.end());
    
    Clear the vectors since they won't be used anymore
    magnetic_field2.clear();
    magnetic_field3.clear();

    std::vector<float> hysteresis_loop(magnetic_field.size());
    int samples = 1000;
    int d = 0;
    for(int sample = 0; sample < samples; ++sample)
    {
        //Do not need to equilibrate for hysteresis loop.  Only perform one iteration of "flips" for each value of h
        initial_lattice = generate_lattice(lattice_size);
        Lattice obj1(lattice_size,T , 1, 0, 1000,10000,initial_lattice);Class initialization
        
        std::vector<float> temp_hysteresis = obj1.get_hysteresis_loop(magnetic_field);
        for(int i = 0; i<magnetic_field.size(); ++i)
        {
            hysteresis_loop[i]+=temp_hysteresis[i]/samples;//Average hysteresis over many samples to get clean data
        }

        std::cout<<float(d)*100/(samples)<<"%\r";
        std::cout.flush();
        d++;
    }
*/


//Uncomment the following to write various data to txt file
/*
    std::ofstream myfile1;
    std::ofstream myfile2;
    std::ofstream myfile3;
    std::ofstream myfile4;
    
    myfile1.open("magnetic_susceptibility.txt");
    myfile2.open("specific_heat.txt");
    myfile3.open("average_magnetization.txt");
    myfile4.open("average_energy.txt");
 for(int size_index = 0; size_index < susceptibility.size(); ++size_index)
    {
        for(int temp_index = 0; temp_index < susceptibility[size_index].size(); ++temp_index)
        {
            myfile1<<susceptibility[size_index][temp_index]<<' ';
            myfile2<<specific_heat[size_index][temp_index]<<' ';
            myfile3<<avg_magnetization[size_index][temp_index]<<' ';
            myfile4<<avg_energy[size_index][temp_index]<<' ';
        }
        myfile1<<"\n";
        myfile2<<"\n";
        myfile3<<"\n";
        myfile4<<"\n";
   }
    myfile1.close();
    myfile2.close();
    myfile3.close();
    myfile4.close();
*/

//Uncomment below to write hysteresis loop data to file
/*
    std::ofstream myfile5;
    myfile5.open("hysteresis_loop.txt");

    for(int h_index = 0; h_index < hysteresis_loop.size(); ++h_index)
    {
        myfile5<<hysteresis_loop[h_index]<<' ';
    }
    myfile5.close();
*/

    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}
