#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdio>

//Domain size
int NX,NY,NUM;

//Time steps
int N=20000;
int NOUTPUT=500;
int NSIGNAL=500;

//Other constants
const int NPOP=9;

double force_x=0.0; 
double force_y=0.0;

//Initial global velocities for hydrodynamics
double global_vel_x=0.01;
double global_vel_y=0.0;

//BGK relaxation parameter

double omega=1.0/0.6;

//Fields and populations
double *f;
double *f2;
double *ux;
double *uy;
double *rho;
double *fx;
double *fy;

//Underlying lattice parameters
double weights[]={4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};
int cx[]={0,1,0,-1,0,1,-1,-1,1};
int cy[]={0,0,1,0,-1,1,1,-1,-1};

//Parameters and structures for particles
struct node_struct {

  /// Constructor

  node_struct() {
   x = 0;
   y = 0;
   x_ref = 0;
   y_ref = 0;
   vel_x = 0;
   vel_y = 0;
   force_x = 0;
   force_y = 0;
  }

  /// Elements

  double x; // current x-position
  double y; // current y-position
  double x_ref; // reference x-position
  double y_ref; // reference y-position
  double vel_x; // node velocity (x-component)
  double vel_y; // node velocity (y-component)
  double force_x; // node force (x-component)
  double force_y; // node force (y-component)
};

struct particle_struct {
    particle_struct() {
        num_nodes = 0;
        radius = 0.0;
        center_x = 0.0;
        center_y = 0.0;
        
        torque = 0.0;
        angle = 0.0;
        omega_rot = 0.0;
        
        force_tot_x = 0.0;
        force_tot_y = 0.0;
        vel_center_x = 0.0;
        vel_center_y = 0.0;
        
        stiffness = 0.03;         
        points = NULL;
        
        rho_ratio = 1.0;
        
        is_moving = false;
     
    }
    
    ~particle_struct(){
        if (num_nodes != 0)
            delete[] points;
    }
    
    int num_nodes;
    double radius, center_x, center_y, vel_center_x, vel_center_y;
    double angle,torque,omega_rot;
    double force_tot_x, force_tot_y;

    double stiffness;
    node_struct* points; 
    
    double rho_ratio;
    
    bool is_moving;       
};

particle_struct* particles;
int num_particles;

void writedensity(std::string const & fname)
{
	std::string filename=fname+".dat";
	std::ofstream fout(filename.c_str());
	fout.precision(10);

	for (int iY=0; iY<NY; iY++)
	{
		for (int iX=0; iX<NX; ++iX)	
		{
			int counter=iY*NX+iX;
			fout<<rho[counter]<<" ";
		}
		fout<<"\n";
	}

}

void writevelocityx(std::string const & fname)
{
	std::string filename=fname+".dat";
	std::ofstream fout(filename.c_str());
	fout.precision(10);

	for (int iY=0; iY<NY; iY++)
	{
		for (int iX=0; iX<NX; ++iX)	
		{
			int counter=iY*NX+iX;
			fout<<ux[counter]<<" ";
		}
		fout<<"\n";
	}

}

void writevelocityy(std::string const & fname)
{
	std::string filename=fname+".dat";
	std::ofstream fout(filename.c_str());
	fout.precision(10);

	for (int iY=0; iY<NY; iY++)
	{
		for (int iX=0; iX<NX; ++iX)	
		{
			int counter=iY*NX+iX;
			fout<<uy[counter]<<" ";
		}
		fout<<"\n";
	}

}

void writevtk(std::string const & fname)
{
	std::string filename=fname+".vtk";
	std::ofstream fout(filename.c_str());
    fout<<"# vtk DataFile Version 3.0\n";
    fout<<"Hydrodynamics representation\n";
    fout<<"ASCII\n\n";
    fout<<"DATASET STRUCTURED_GRID\n";
    fout<<"DIMENSIONS "<<NX<<" "<<NY<<" "<<1<<"\n";
    fout<<"POINTS "<<NX*NY*1<<" double\n";
    for(int counter=0;counter<NUM;counter++)
    {
    	int iX=counter%NX;
    	int iY=counter/NX;
        fout<<iX<<" "<<iY-0.5<<" "<<0<<"\n";
    }
    fout<<"\n";
    fout<<"POINT_DATA "<<NX*NY*1<<"\n";
    
    fout<<"SCALARS density double\n";
    fout<<"LOOKUP_TABLE density_table\n";
    for(int counter=0;counter<NUM;counter++)
   		fout<<rho[counter]<<"\n";
   		
    fout<<"VECTORS velocity double\n";
    for(int counter=0;counter<NUM;counter++)
    		fout<<ux[counter]<<" "<<uy[counter]<<" 0\n";
    
    // Now write all the immersed boundary nodes
    int num_overall = 0;
    for (int iPart = 0; iPart < num_particles; iPart++)
    {
        num_overall = num_overall + particles[iPart].num_nodes;
    }
    
    std::string filename_imm=fname+"_imm.vtk";
    std::ofstream fout_imm(filename_imm.c_str());		
    //Immersed boundary datasets
    fout_imm << "# vtk DataFile Version 3.0 \n";
    fout_imm << "Immersed boundary representation\n";
    fout_imm << "ASCII\n\n";
    fout_imm << "DATASET POLYDATA\n";
  	fout_imm << "POINTS " << num_overall << " double\n";

    for (int iPart = 0; iPart < num_particles; iPart++)
  	    for(int n = 0; n < particles[iPart].num_nodes; n++) 
    	    fout_imm << particles[iPart].points[n].x << " " << particles[iPart].points[n].y << " 0\n";
  	

  	//Write lines between neighboring nodes
  	fout_imm << "LINES " << num_overall << " " << 3 * num_overall << "\n";
    
    int num_curr = 0;
    for (int iPart = 0; iPart < num_particles; iPart++)
    {
        int num_curr_beg = num_curr;
  	    for(int n = 0; n < particles[iPart].num_nodes - 1; n++)
  	    {    
    	    fout_imm << "2 " << num_curr << " " << num_curr + 1 << "\n";
    	    num_curr++;
        }
        fout_imm << "2 " << num_curr << " " << num_curr_beg << "\n";
        num_curr++;
     }   
	
	//Write vertices
  	fout_imm << "VERTICES 1 " << num_overall + 1 << "\n";
  	fout_imm << num_overall <<" ";

  	for(int n = 0; n < num_overall; n++)
    	fout_imm << n << " ";
    fout_imm.close();
    

    std::string filename_ref=fname+"_ref.vtk";
    std::ofstream fout_ref(filename_ref.c_str());		
    //Immersed boundary datasets
    fout_ref << "# vtk DataFile Version 3.0 \n";
    fout_ref << "Immersed boundary representation\n";
    fout_ref << "ASCII\n\n";
    fout_ref << "DATASET POLYDATA\n";
  	fout_ref << "POINTS " << num_overall << " double\n";

    for (int iPart = 0; iPart < num_particles; iPart++)
  	    for(int n = 0; n < particles[iPart].num_nodes; n++) 
    	    fout_ref << particles[iPart].points[n].x_ref << " " << particles[iPart].points[n].y_ref << " 0\n";
  	

  	//Write lines between neighboring nodes
  	fout_ref << "LINES " << num_overall << " " << 3 * num_overall << "\n";
    
    num_curr = 0;
    for (int iPart = 0; iPart < num_particles; iPart++)
    {
        int num_curr_beg = num_curr;
  	    for(int n = 0; n < particles[iPart].num_nodes - 1; n++)
  	    {    
    	    fout_ref << "2 " << num_curr << " " << num_curr + 1 << "\n";
    	    num_curr++;
        }
        fout_ref << "2 " << num_curr << " " << num_curr_beg << "\n";
        num_curr++;
     }   
	
	//Write vertices
  	fout_ref << "VERTICES 1 " << num_overall + 1 << "\n";
  	fout_ref << num_overall <<" ";

  	for(int n = 0; n < num_overall; n++)
    	fout_ref << n << " ";
    fout_ref.close();  

        
}

void init_hydro()
{
    NY = 82;
    NX = 400;
    NUM=NX*NY;
    rho=new double[NUM];
    ux=new double[NUM];
    uy=new double[NUM];
    fx=new double[NUM];
    fy=new double[NUM];
 
 	//Initialization of density
    for(int counter=0;counter<NUM;counter++)
    {
		int iY = counter / NX;			
		ux[counter]=global_vel_x; 
		uy[counter]=global_vel_y;
		fx[counter]=0.0;
		fy[counter]=0.0;
		rho[counter]=1.0;
	}
	for(int counter=0; counter<NUM; counter++)
	{
		int iY = counter / NX;
		int iX = counter % NX;
	    
	    for(int iPart = 0; iPart < num_particles; iPart++)
	    {
	        double cx = particles[iPart].center_x;
	        double cy = particles[iPart].center_y;
	        double rad = particles[iPart].radius;
	        	
		    if ((iY-cy)*(iY-cy)+(iX-cx)*(iX-cx) <= rad*rad)
		    {
			    ux[counter]=0.0;
			    uy[counter]=0.0;
		    }
		}
	}
}

void init_immersed()
{
    // Here we initialize all the particles
    num_particles = 2; 
    particles = new particle_struct[num_particles];
	
    particles[0].radius = 20;
    particles[0].center_x = 50;
    particles[0].center_y = 30;
    particles[0].is_moving = true;
    
    particles[1].radius = 20;
    particles[1].center_x = 200;
    particles[1].center_y = 40;
    particles[1].is_moving = false;    
    
    for(int iPart = 0; iPart < num_particles; iPart++)
    {
        particles[iPart].num_nodes = 150;
        int num = particles[iPart].num_nodes;        
        particles[iPart].points = new node_struct[num];
        
        double cx = particles[iPart].center_x;
        double cy = particles[iPart].center_y;
        double rad = particles[iPart].radius;
        for(int n = 0; n < num; ++n) 
        {
            particles[iPart].points[n].x = cx + rad * cos(2.0 * M_PI * double(n) / double(num));
            particles[iPart].points[n].x_ref = cx + rad * cos(2.0 * M_PI * double(n) / double(num));
            particles[iPart].points[n].y = cy + rad * sin(2.0 * M_PI * double(n) / double(num));
            particles[iPart].points[n].y_ref = cy + rad * sin(2.0 * M_PI * double(n) / double(num));
        }
    }
   
   
   
}

void init()
{
	//Creating arrays
	f =new double[NUM*NPOP];
	f2=new double[NUM*NPOP]; 
    
	//Bulk nodes initialization
	double feq;
	
	for(int iY=0;iY<NY;iY++)
		for(int iX=0;iX<NX;iX++)
		{
			int  counter=iY*NX+iX;
			double dense_temp=rho[counter];
			double ux_temp=ux[counter];
			double uy_temp=uy[counter];
            for (int k=0; k<NPOP; k++)
			{
				feq=weights[k]*(dense_temp+3.0*dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp
								                +(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp
								                +2.0*ux_temp*uy_temp*cx[k]*cy[k]));
                f[counter*NPOP+k]=feq;
			}
			
		}

}

void compute_particle_forces()
{
  //Reset forces
    for(int iPart = 0; iPart < num_particles; iPart++)
    {
        int num = particles[iPart].num_nodes;
        for(int n = 0; n < num; n++) 
        {   
            particles[iPart].points[n].force_x = 0.0;
            particles[iPart].points[n].force_y = 0.0;
        }
        //Compute rigid forces
        double area = 2.0 * M_PI * particles[iPart].radius / num; 
        for(int n = 0; n < num; n++) 
        {
  	        // Don't understand why we need area?
            particles[iPart].points[n].force_x = -particles[iPart].stiffness * 
            (particles[iPart].points[n].x - particles[iPart].points[n].x_ref) * area;
            
            particles[iPart].points[n].force_y = -particles[iPart].stiffness * 
            (particles[iPart].points[n].y - particles[iPart].points[n].y_ref) * area;
        }
  
        particles[iPart].torque = 0.0;
        particles[iPart].force_tot_x = 0.0;
        particles[iPart].force_tot_y = 0.0;
        for(int n = 0; n < num; n++)
        {
  	        double fx = particles[iPart].points[n].force_x;
  	        double fy = particles[iPart].points[n].force_y;
  	        double rx = particles[iPart].points[n].x - particles[iPart].center_x;
  	        double ry = particles[iPart].points[n].y - particles[iPart].center_y;
  	        particles[iPart].torque += fx*ry - fy*rx;
  	        particles[iPart].force_tot_x += fx;
  	        particles[iPart].force_tot_y += fy;
        }  	  
    }
  
}

void spread_particle_forces()
{
    // Reset forces
    for(int iX = 0; iX < NX; iX++) {
        for(int iY = 0; iY < NY; iY++) {
            int counter=iY*NX+iX;
            fx[counter] = 0;
            fy[counter] = 0;
        }
    }

    for (int iPart = 0; iPart < num_particles; iPart++)
    {
        int num = particles[iPart].num_nodes;
        // Spread forces
        for(int n = 0; n < num; n++) 
        {

            // Identify the lowest fluid lattice node in interpolation range.
            // 'Lowest' means: its x- and y-values are the smallest.
            // The other fluid nodes in range have coordinates
            // (x_int + 1, y_int), (x_int, y_int + 1), and (x_int + 1, y_int + 1).
            double x = particles[iPart].points[n].x;
            double y = particles[iPart].points[n].y;
            
            int xset[2], xset2[2];
	        int yset[2], yset2[2];
	        
	        xset[0] = (int(floor(x)) % NX + NX) % NX;
	        xset[1] = (int(ceil(x)) % NX + NX) % NX;
	        yset[0] = (int(floor(y)) % NY + NY) % NY;
	        yset[1] = (int(ceil(y)) % NY + NY) % NY;
	        
	        xset2[0] = int(floor(x));
	        xset2[1] = int(ceil(x));
	        yset2[0] = int(floor(y));
	        yset2[1] = int(ceil(y));    
	        
            for(int i = 0; i < 2; i++) 
            {
    	        for(int j = 0; j < 2; j++) 
      	        {

                	// Compute distance between object node and fluid lattice node.
    
 	               	double dist_x = fabs(x-xset2[i]);
        		    double dist_y = fabs(y-yset2[j]);

                	// Compute interpolation weights for x- and y-direction based on the distance.

                	double weight_x = 1.0 - dist_x;
        	        double weight_y = 1.0 - dist_y;
        	        int counter = yset[j] * NX + xset[i];

        	        // Compute lattice force.
			        fx[counter] += particles[iPart].points[n].force_x * weight_x * weight_y;
			        fy[counter] += particles[iPart].points[n].force_y * weight_x * weight_y;

                }
            }
        }
    }
        
}

void interpolate_particle_velocities()
{


    for (int iPart = 0; iPart < num_particles; iPart++)
    {
        int num = particles[iPart].num_nodes;
        
        for(int n = 0; n < num; n++) 
        {

            particles[iPart].points[n].vel_x = 0;
            particles[iPart].points[n].vel_y = 0;
            
            double x = particles[iPart].points[n].x;
            double y = particles[iPart].points[n].y;

            // Identify the lowest fluid lattice node in interpolation range (see spreading).
            int xset[2], xset2[2];
	        int yset[2], yset2[2];
	        
	        xset[0] = (int(floor(x)) % NX + NX) % NX;
	        xset[1] = (int(ceil(x)) % NX + NX) % NX;
	        yset[0] = (int(floor(y)) % NY + NY) % NY;
	        yset[1] = (int(ceil(y)) % NY + NY) % NY;
	        
	        xset2[0] = int(floor(x));
	        xset2[1] = int(ceil(x));
	        yset2[0] = int(floor(y));
	        yset2[1] = int(ceil(y));
            
            for(int i = 0; i < 2; i++) 
            {
                for(int j = 0; j < 2; j++) 
                {

                    // Compute distance between object node and fluid lattice node.

                    double dist_x = fabs(x - xset2[i]);
                    double dist_y = fabs(y - yset2[j]);

                    // Compute interpolation weights for x- and y-direction based on the distance.
                    double weight_x = 1.0 - dist_x;
                    double weight_y = 1.0 - dist_y;

                    // Compute node velocities.
		            int counter =yset[j] * NX + xset[i];
		
                    particles[iPart].points[n].vel_x += ux[counter] * weight_x * weight_y;
                    particles[iPart].points[n].vel_y += uy[counter] * weight_x * weight_y;
                }
            }
        }
  
    }

}

void update_particle_position() 
{
    for(int iPart=0; iPart < num_particles; iPart++)
    {
    
        //Update node and center positions  
        int num = particles[iPart].num_nodes;
        double rad = particles[iPart].radius;
        double ratio = particles[iPart].rho_ratio;
        double torque = particles[iPart].torque;
        
        for(int n = 0; n < num; n++) 
        {
            double x = particles[iPart].points[n].x;
            double y = particles[iPart].points[n].y;
            x = x + particles[iPart].points[n].vel_x;
            y = y + particles[iPart].points[n].vel_y;
            particles[iPart].points[n].x = x;
            particles[iPart].points[n].y = y;
        }

        // Account for motion 
        
        if (particles[iPart].is_moving)
        {
            particles[iPart].vel_center_x = particles[iPart].vel_center_x 
                                          - particles[iPart].force_tot_x/(M_PI*rad*rad*ratio);
            particles[iPart].vel_center_y = particles[iPart].vel_center_y 
                                          - particles[iPart].force_tot_y/(M_PI*rad*rad*ratio);
                                      
            particles[iPart].center_x = particles[iPart].center_x 
                                      + particles[iPart].vel_center_x;
            particles[iPart].center_y = particles[iPart].center_y 
                                      + particles[iPart].vel_center_y;
  
            particles[iPart].omega_rot = particles[iPart].omega_rot 
                                       - 2.0*torque/(ratio*M_PI*rad*rad*rad*rad);
            particles[iPart].angle = particles[iPart].angle + particles[iPart].omega_rot;

            for(int n = 0; n < num; ++n) 
            {
                double x = rad*cos(2.0 * M_PI * double(n) / double(num));
                double y = rad*sin(2.0 * M_PI * double(n) / double(num));
                double angle = particles[iPart].angle;
                particles[iPart].points[n].x_ref = particles[iPart].center_x + x*cos(angle) + y*sin(angle);
                particles[iPart].points[n].y_ref = particles[iPart].center_y - x*sin(angle) + y*cos(angle);
            }
        }
    
    }

    return;
}

void collide_bulk()
{

    for(int counter=0;counter<NUM;counter++)
    {			
		int iX=counter%NX;
		int iY=counter/NX;
		
		double dense_temp=0.0;
		double ux_temp=0.0;
		double uy_temp=0.0;
		       
		for(int k=0;k<NPOP;k++)
		{
			dense_temp+=f[counter*NPOP+k];
			ux_temp+=f[counter*NPOP+k]*cx[k];
			uy_temp+=f[counter*NPOP+k]*cy[k];
		}
		// Here needs to be a reference of forces
		ux_temp=(ux_temp+0.5*(fx[counter]+force_x))/dense_temp;
		uy_temp=(uy_temp+0.5*(fy[counter]+force_y))/dense_temp;
        	
		rho[counter]=dense_temp;
		ux[counter]=ux_temp;
		uy[counter]=uy_temp;
	        
		double sum=0.0;
        
      	double feqeq[NPOP],force[NPOP];
      	
		for (int k=1; k<NPOP; k++)
		{
			feqeq[k]=weights[k]*(dense_temp+3.0*dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
                	+4.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp
                	+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+2.0*ux_temp*uy_temp*cx[k]*cy[k]));
            sum+=feqeq[k];
		}

		feqeq[0]=dense_temp-sum;
        
 		//Obtain force population
        for (int k=0;k<NPOP;k++)
        {
			force[k]=weights[k]*(1.0-0.5*omega)*
					 (3.0*(fx[counter]+force_x)*cx[k]+3.0*(fy[counter]+force_y)*cy[k]+
        	         9.0*((cx[k]*cx[k]-1.0/3.0)*(fx[counter]+force_x)*ux_temp+
        	               cx[k]*cy[k]*((fx[counter]+force_x)*uy_temp+(fy[counter]+force_y)*ux_temp)+
 						  (cy[k]*cy[k]-1.0/3.0)*(fy[counter]+force_y)*uy_temp));
        }
		
		for(int k=0; k < NPOP; k++)
		{
			f2[counter*NPOP+k]=f[counter*NPOP+k]*(1.0-omega)+omega*feqeq[k]+force[k];    	
        }
    }

}

void update_inamuro()
{
    //Perform velocity bounce back on the inlet and outlet
    for(int iY=0;iY<NY;iY++)
    {
    	int offset = iY*NX*NPOP;
    	double rho_wall  = 1.0/(1.0-global_vel_x)*(f2[offset+0]+f2[offset+2]+f2[offset+4]+
    	                   2.0*(f2[offset+3]+f2[offset+6]+f2[offset+7]));
        double rho_prime = 6.0/(1.0+3.0*global_vel_x+3.0*global_vel_x*global_vel_x)*
                           (rho_wall*global_vel_x+(f2[offset+3]+f2[offset+6]+f2[offset+7]));
    	
    	f2[offset+1]=weights[1]*rho_prime*(1.0+3.0*global_vel_x+3.0*global_vel_x*global_vel_x);
    	f2[offset+5]=weights[5]*rho_prime*(1.0+3.0*global_vel_x+3.0*global_vel_x*global_vel_x);
    	f2[offset+8]=weights[8]*rho_prime*(1.0+3.0*global_vel_x+3.0*global_vel_x*global_vel_x);

        // Recalculate macroscopic parameters
        rho[iY*NX] = 0.0;
        ux[iY*NX] = 0.0;
        uy[iY*NX] = 0.0;
        for (int iPop=0;iPop<NPOP;iPop++)
        {
            rho[iY*NX] += f2[offset+iPop];
            ux[iY*NX] += f2[offset+iPop]*cx[iPop];
            uy[iY*NX] += f2[offset+iPop]*cy[iPop];
        }
        ux[iY*NX]=ux[iY*NX]/rho[iY*NX];
        uy[iY*NX]=uy[iY*NX]/rho[iY*NX];

        offset = (iY*NX+NX-1)*NPOP;
    	rho_wall  = 1.0/(1.0+global_vel_x)*(f2[offset+0]+f2[offset+2]+f2[offset+4]+
                    2.0*(f2[offset+1]+f2[offset+5]+f2[offset+8]));
        rho_prime = 6.0/(1.0-3.0*global_vel_x+3.0*global_vel_x*global_vel_x)*
                    (-rho_wall*global_vel_x+(f2[offset+1]+f2[offset+5]+f2[offset+8]));
    	
    	f2[offset+3]=weights[3]*rho_prime*(1.0-3.0*global_vel_x+3.0*global_vel_x*global_vel_x);
    	f2[offset+6]=weights[6]*rho_prime*(1.0-3.0*global_vel_x+3.0*global_vel_x*global_vel_x);
    	f2[offset+7]=weights[7]*rho_prime*(1.0-3.0*global_vel_x+3.0*global_vel_x*global_vel_x);

        // Recalculate macroscopic parameters
        rho[iY*NX+NX-1] = 0.0;
        ux[iY*NX+NX-1] = 0.0;
        uy[iY*NX+NX-1] = 0.0;
        for (int iPop=0;iPop<NPOP;iPop++)
        {
            rho[iY*NX+NX-1] += f2[offset+iPop];
            ux[iY*NX+NX-1] += f2[offset+iPop]*cx[iPop];
            uy[iY*NX+NX-1] += f2[offset+iPop]*cy[iPop];
        }
        
        ux[iY*NX+NX-1]=ux[iY*NX+NX-1]/rho[iY*NX+NX-1];
        uy[iY*NX+NX-1]=uy[iY*NX+NX-1]/rho[iY*NX+NX-1];
    	
    }
    
}

// The code below is not used
#if 0
void update_bounce_back()
{
    //Perform velocity bounce back on both walls
    for(int iX=0;iX<NX;iX++)
    {
        int iX_top    = (iX + 1 + NX ) % NX;
        int iX_bottom = (iX - 1 + NX ) % NX;
    	f2[((NY-1)*NX+iX)*NPOP+4]=f2[((NY-2)*NX+iX)*NPOP+2]       +2.0*weights[4]*3.0*(-vel_wall*cx[4]);
        f2[((NY-1)*NX+iX)*NPOP+7]=f2[((NY-2)*NX+iX_bottom)*NPOP+5]+2.0*weights[7]*3.0*(-vel_wall*cx[7]); 
    	f2[((NY-1)*NX+iX)*NPOP+8]=f2[((NY-2)*NX+iX_top)*NPOP+6]   +2.0*weights[8]*3.0*(-vel_wall*cx[8]);

    	f2[iX*NPOP+2]=f2[(NX+iX)*NPOP+4]       +2.0*weights[2]*3.0*(vel_wall*cx[2]);
        f2[iX*NPOP+5]=f2[(NX+iX_top)*NPOP+7]   +2.0*weights[5]*3.0*(vel_wall*cx[5]); 
    	f2[iX*NPOP+6]=f2[(NX+iX_bottom)*NPOP+8]+2.0*weights[6]*3.0*(vel_wall*cx[6]);
    }

}
#endif

void finish_simulation()
{
	delete[] rho;
	delete[] ux;
	delete[] uy;
	delete[] f;
	delete[] f2;
}

void stream()
{
    for(int counter=0;counter<NUM;counter++)
	{
		int iX=counter%NX;
		int iY=counter/NX;

		for(int iPop=0;iPop<NPOP;iPop++)
		{
			int iX2=(iX-cx[iPop]+NX)%NX;
			int iY2=(iY-cy[iPop]+NY)%NY;
			int counter2=iY2*NX+iX2;
			f[counter*NPOP+iPop]=f2[counter2*NPOP+iPop];
		}
	}

	
}

void calculate_overall_force(std::string filename,int counter)
{

    std::cout << "***** New signal time *****" << std::endl;
    std::cout << "Time: " << counter << std::endl;
    for (int iPart = 0; iPart < num_particles; iPart++)
    {
        std::cout << "Particle: "  << iPart << std::endl;
        double rad = particles[iPart].radius;
        int num = particles[iPart].num_nodes;
        double len = 2.0*M_PI*rad/double(num);
        std::cout << "Torque: " << particles[iPart].torque << std::endl;
        std::cout << "Angle: " << particles[iPart].angle << std::endl;
  	    std::cout << "Center_x: " << particles[iPart].center_x << std::endl;
  	    std::cout << "Center_y: " << particles[iPart].center_y << std::endl;
        std::cout << "Vel_x:" << particles[iPart].vel_center_x << std::endl;        
        std::cout << "Vel_y:" << particles[iPart].vel_center_y << std::endl;  

        double aver_deformation = 0.0;
        for(int n = 0; n < num; n++) 
        { 
            double x     = particles[iPart].points[n].x;
            double x_ref = particles[iPart].points[n].x_ref;
            double y     = particles[iPart].points[n].y;
            double y_ref = particles[iPart].points[n].y_ref;
            aver_deformation += sqrt((x-x_ref)*(x-x_ref)+(y-y_ref)*(y-y_ref));
        }
        aver_deformation = aver_deformation/num;
       
        std::cout << "Len: " << len << std::endl;
        std::cout << "Average deformation: " << aver_deformation << std::endl;

    }
    
}

int main(int argc, char* argv[])
{
	init_immersed();	
	init_hydro();
    init();

	for(int counter=0;counter<=N;counter++)
	{
		compute_particle_forces();
		spread_particle_forces();
        collide_bulk();
        update_inamuro();
		stream();
		interpolate_particle_velocities();
		update_particle_position();
        
	    if (counter%NSIGNAL==0)
	    {
	    	std::cout<<"Counter="<<counter<<"\n";
	    	std::stringstream filewriteoverall;
	    	filewriteoverall << "coors.dat";
			calculate_overall_force(filewriteoverall.str(),counter);	    	
    	}
    	
		//Writing files
		if (counter%NOUTPUT==0)
		{
            std::cout<<"Output "<<counter<<"\n";
			std::stringstream filewritedensity;
  			std::stringstream filewritevelocityx;
  			std::stringstream filewritevelocityy;
 			std::stringstream filevtk;
 			
 			std::stringstream counterconvert;
 			counterconvert<<counter;
 			filewritedensity<<std::fixed;
 			filewritevelocityx<<std::fixed;
 			filewritevelocityy<<std::fixed;
 			filevtk<<std::fixed;

			//filewritedensity<<"den"<<std::string(7-counterconvert.str().size(),'0')<<counter;
			filewritevelocityx<<"velx"<<std::string(7-counterconvert.str().size(),'0')<<counter;
			//filewritevelocityy<<"vely"<<std::string(7-counterconvert.str().size(),'0')<<counter;
			filevtk<<"vtk"<<std::string(7-counterconvert.str().size(),'0')<<counter;
			
 			//writedensity(filewritedensity.str());
 			writevelocityx(filewritevelocityx.str());
 			//writevelocityy(filewritevelocityy.str());
 			writevtk(filevtk.str());
		}

	}
    
   	return 0;
}
