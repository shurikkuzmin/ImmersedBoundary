#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

//Domain size
int NX,NY,NUM;

//Time steps
int N=1000;
int NOUTPUT=100;
int NSIGNAL=100;

//Other constants
const int NPOP=9;

const double rho_wall=0.5;

double force_x=0.000006;
double force_y=0.0;

//BGK relaxation parameter
double omega=1.0; // /0.525;


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

//IB particle values
int NUMIB;
node_struct* points;
double center_x = 30.0;
double center_y = 30.0;
double radius   = 10.0;

//Global stiffness
double stiffness = 0.003;

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
        fout<<iX<<" "<<iY<<" "<<0<<"\n";
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
    
    std::string filename_imm=fname+"_imm.vtk";
    std::ofstream fout_imm(filename_imm.c_str());		
    //Immersed boundary datasets
    fout_imm << "# vtk DataFile Version 3.0 \n";
    fout_imm << "Immersed boundary representation\n";
    fout_imm << "ASCII\n\n";
    fout_imm << "DATASET POLYDATA\n";
  	fout_imm << "POINTS " << NUMIB << " double\n";

  	for(int n = 0; n < NUMIB; n++) 
    	fout_imm << points[n].x << " " << points[n].y << " 0\n";
  	

  	//Write lines between neighboring nodes
  	fout_imm << "LINES " << NUMIB << " " << 3 * NUMIB << "\n";

  	for(int n = 0; n < NUMIB; n++) 
    	fout_imm << "2 " << n << " " << (n + 1) % NUMIB << "\n";

	//Write vertices

  	fout_imm << "VERTICES 1 " << NUMIB + 1 << "\n";
  	fout_imm << NUMIB << " ";

  	for(int n = 0; n < NUMIB; n++) 
    	fout_imm << n << " ";
  
    fout_imm.close();
        
}

void init_hydro()
{
    NY=62;
    NX=220;
    NUM=NX*NY;
    rho=new double[NUM];
    ux=new double[NUM];
    uy=new double[NUM];
    fx=new double[NUM];
    fy=new double[NUM];
 
 	//Initialization of density
    for(int counter=0;counter<NUM;counter++)
    {
		ux[counter]=0.0;
		uy[counter]=0.0;
		fx[counter]=0.0;
		fy[counter]=0.0;
		rho[counter]=1.0;
	}
    //Initialization of density walls
    for(int iX=0;iX<NX;iX++)
    {
        rho[iX]=rho_wall;
        rho[NX*(NY-1)+iX]=rho_wall;
    }
    
}

void init_immersed()
{
	//Number of immersed boundary points. We take it as perimeter with 
	NUMIB = 80;
	points = new node_struct[NUMIB];	

    for(int n = 0; n < NUMIB; ++n) {
        points[n].x = center_x + radius * cos(2.0 * M_PI * double(n) / double(NUMIB));
        points[n].x_ref = center_x + radius * cos(2.0 * M_PI * double(n) / double(NUMIB));
        points[n].y = center_y + radius * sin(2.0 * M_PI * double(n) / double(NUMIB));
        points[n].y_ref = center_y + radius * sin(2.0 * M_PI * double(n) / double(NUMIB));
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
  for(int n = 0; n < NUMIB; n++) {
    points[n].force_x = 0.0;
    points[n].force_y = 0.0;
  }
  //Compute rigid forces
  double area = 2.0 * M_PI * radius / NUMIB; 
  for(int n = 0; n < NUMIB; n++) 
  {
      points[n].force_x = -stiffness * (points[n].x - points[n].x_ref) * area;
      points[n].force_y = -stiffness * (points[n].y - points[n].y_ref) * area;
  }
	
}

void spread_particle_forces()
{
  // Reset forces
  
  for(int iX = 0; iX < NX; iX++) {
    for(int iY = 1; iY < NY - 1; iY++) {
      int counter=iY*NX+iX;
      fx[counter] = 0;
      fy[counter] = 0;
    }
  }

  // Spread forces
  
  for(int n = 0; n < NUMIB; n++) {

    // Identify the lowest fluid lattice node in interpolation range.
    // 'Lowest' means: its x- and y-values are the smallest.
    // The other fluid nodes in range have coordinates
    // (x_int + 1, y_int), (x_int, y_int + 1), and (x_int + 1, y_int + 1).

	int xbottom = floor(points[n].x);
	int ybottom = floor(points[n].y);
	int xtop    = ceil(points[n].x); 
	int ytop    = ceil(points[n].y);
	std::cout << "X: " << points[n].x << " " << xbottom << " " << xtop << "\n";
	std::cout << "Y: " << points[n].y << " " << ybottom << " " << ytop << "\n";
	if(xbottom == xtop)
	{
		std::cout << "Zhopa  xbottom=xtop" << n << "\n"; 
	}
	if(ybottom == ytop)
	{
		std::cout << "Zhopa  ybottom=ytop" << n << "\n"; 
	}
 
    for(int iX = xbottom; iX <= xtop; iX++) 
    {
    	for(int iY = ybottom; iY <= ytop; iY++) 
      	{

        	// Compute distance between object node and fluid lattice node.

 	       	double dist_x = abs(points[n].x-iX);
    		double dist_y = abs(points[n].y-iY);

        	// Compute interpolation weights for x- and y-direction based on the distance.

        	double weight_x = 1.0 - dist_x;
        	double weight_y = 1.0 - dist_y;
        	
        	int counter = iX * NX + iY;

        	// Compute lattice force.
			fx[counter] += points[n].force_x * weight_x * weight_y;
			fy[counter] += points[n].force_y * weight_x * weight_y; 
      }
    }
  }
}

void interpolate_particle_velocities()
{

  for(int n = 0; n < NUMIB; n++) {

    points[n].vel_x = 0;
    points[n].vel_y = 0;

    // Identify the lowest fluid lattice node in interpolation range (see spreading).

    //int x_int = (int) (particle.node[n].x - 0.5 + Nx) - Nx;
    //int y_int = (int) (particle.node[n].y + 0.5);
    
    int xbottom = floor(points[n].x);
    int ybottom = floor(points[n].y);
    int xtop    = ceil(points[n].x);
    int ytop    = ceil(points[n].y);
    
    for(int iX = xbottom; iX <= xtop; iX++) 
    {
      for(int iY = ybottom; iY <= ytop; iY++) 
      {

        // Compute distance between object node and fluid lattice node.

        double dist_x = abs(points[n].x - iX);
        double dist_y = abs(points[n].y - iY);

        // Compute interpolation weights for x- and y-direction based on the distance.

        double weight_x = 1 - dist_x;
        double weight_y = 1 - dist_y;

        // Compute node velocities.
		int counter = iY * NX + iX;
		
        points[n].vel_x += ux[counter] * weight_x * weight_y;
        points[n].vel_y += uy[counter] * weight_x * weight_y;
      }
    }
  }

  return;
}

void update_particle_position() 
{
  //We don't actually update the center of the particle as it is supposed
  //to be fixed
  
  //Reset center position
  //center_x = 0;
  //center_y = 0;

  //Update node and center positions

  for(int n = 0; n < NUMIB; n++) {
    points[n].x = fmod(points[n].x + points[n].vel_x + double(NX),double(NX));
    points[n].y = points[n].y + points[n].vel_y;
    //particle.center.x += particle.node[n].x / particle.num_nodes;
    //particle.center.y += particle.node[n].y / particle.num_nodes;
  }

  return;
}

void collide_bulk()
{

    for(int counter=0;counter<NUM;counter++)
    {			
		int iX=counter%NX;
		int iY=counter/NX;
		
		if ( (iY==0) || (iY==NY-1) )
            continue;
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

void update_bounce_back()
{
	//Perform bounce back on walls
    for(int iX=0;iX<NX;iX++)
    {
        int iX_top=(iX+1+NX)%NX;
        int iX_bottom=(iX-1+NX)%NX;
        
        f2[iX*NPOP+2]=f2[(NX+iX)*NPOP+4];
        f2[iX*NPOP+5]=f2[(NX+iX_top)*NPOP+7];
        f2[iX*NPOP+6]=f2[(NX+iX_bottom)*NPOP+8];
        f2[(NX*(NY-1)+iX)*NPOP+4]=f2[(NX*(NY-2)+iX)*NPOP+2];
        f2[(NX*(NY-1)+iX)*NPOP+7]=f2[(NX*(NY-2)+iX_bottom)*NPOP+5];
        f2[(NX*(NY-1)+iX)*NPOP+8]=f2[(NX*(NY-2)+iX_top)*NPOP+6];        
    }
}

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

        if ((iY==0) || (iY==NY-1))
            continue;
		for(int iPop=0;iPop<NPOP;iPop++)
		{
			int iX2=(iX-cx[iPop]+NX)%NX;
			int iY2=(iY-cy[iPop]+NY)%NY;
			int counter2=iY2*NX+iX2;
			f[counter*NPOP+iPop]=f2[counter2*NPOP+iPop];
		}
	}

	
}

void calculate_overall_force()
{

  double force_tot_x = 0.0;
  double force_tot_y = 0.0;
  double vel_center_x = 0.0;
  double vel_center_y = 0.0;
  double aver_deformation = 0.0;

  for(int n = 0; n < NUMIB; n++) 
  {
    force_tot_x  += points[n].force_x;
    force_tot_y  += points[n].force_y;
    vel_center_x += points[n].vel_x;
    vel_center_y += points[n].vel_y;
    aver_deformation += sqrt((points[n].x-points[n].x_ref)*(points[n].x-points[n].x_ref)+
    						 (points[n].y-points[n].y_ref)*(points[n].y-points[n].y_ref));
  }
  
  aver_deformation = aver_deformation/NUMIB;
  
  std::cout << "Drag force: " << force_tot_x << "\n"; 
  std::cout << "Lift force: " << force_tot_y << "\n"; 
  std::cout << "Vel center X: " << vel_center_x << "\n";
  std::cout << "Vel center Y: " << vel_center_y << "\n";
  std::cout << "Average deformation: " << aver_deformation << "\n";
  std::cout << "Average scaled force: " << aver_deformation*NUMIB << "\n";
}

int main(int argc, char* argv[])
{

	
	init_hydro();
	init_immersed();
    init();

	for(int counter=0;counter<=N;counter++)
	{
		std::cout << "Counter=" << counter << "\n";
		compute_particle_forces();
		spread_particle_forces();
        collide_bulk();
        update_bounce_back();
		stream();
		interpolate_particle_velocities();
		update_particle_position();
        
	    if (counter%NSIGNAL==0)
	    {
	    	std::cout<<"Counter="<<counter<<"\n";
    	}
    	
		//Writing files
		if (counter%NOUTPUT==0)
		{
			calculate_overall_force();
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
