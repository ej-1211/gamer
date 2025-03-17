#include "GAMER.h"
#include "TestProb.h"
#include "Par_EquilibriumIC.h"
#include "string"
#include <stdio.h>
#include "global_var.h"
#include <cmath>
using namespace std;

// negligibly small uniform density and energy
double GC_SmallGas;

// parameter for center setting
bool FixCenter;
double SearchRadius;

double GC_rr;
double GC_theta;
double GC_mm;
// GC position
double GC_xx = 100000000.0;
double GC_yy = 100000000.0;
double GC_zz = 100000000.0;

// Initial Halo center
double Halo_Center_x_Initial;
double Halo_Center_y_Initial;
double Halo_Center_z_Initial;

// for RESTART
double previous_center_x;
double previous_center_y;
double previous_center_z;

// declare the potential minimum last step
double min_pot_last[3] ;

// declare the output profile's dt and the ideal output t
double ideal_output_t = 0.0;
double output_profile_dt;

// declare the radius and density information
// declare a global variable for the profile
Profile_t *Prof_[1];


//// problem-specific function prototypes


void Mis_UserWorkBeforeNextLevel_Find_GC( const int lv, const double TimeNew, const double TimeOld, const double dt );

void AdjustGCPotential(const bool add, const double GC_x, const double GC_y, const double GC_z);
void User_Output();
// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Aux_Record_User_GC();
// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Init_User_moveGC();


//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


// errors
#  if ( MODEL != ELBDM )
   Aux_Error( ERROR_INFO, "MODEL != ELBDM !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
//   if ( OPT__INIT == INIT_BY_FUNCTION  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
//     Aux_Error( ERROR_INFO, "please set PAR_INIT = 1 (by FUNCTION) !!\n" );
#  else
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n" );
#  endif

#  ifndef SUPPORT_GSL
   Aux_Error( ERROR_INFO, "SUPPORT_GSL must be enabled !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
      for (int f=0; f<6; f++)
      if ( OPT__BC_FLU[f] == BC_FLU_PERIODIC )
         Aux_Message( stderr, "WARNING : periodic BC for fluid is not recommended for this test !!\n" );

#     ifdef GRAVITY
      if ( OPT__BC_POT == BC_POT_PERIODIC )
         Aux_Message( stderr, "WARNING : periodic BC for gravity is not recommended for this test !!\n" );
#     endif
   } // if ( MPI_Rank == 0 )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == ELBDM && defined GRAVITY )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    =  10;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" ); }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }

// load run-time parameters
   const char* FileName = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;
   // ********************************************************************************************************************************
   // ReadPara->Add( "KEY_IN_THE_FILE",      &VARIABLE,              DEFAULT,       MIN,              MAX               );
   // ********************************************************************************************************************************
   ReadPara->Add( "GC_SmallGas",             &GC_SmallGas,           1e-10,          0.,               NoMax_double      );

   ReadPara->Add( "GC_POSX",                 &GC_xx,               NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "GC_POSY",                 &GC_yy,               NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "GC_POSZ",                 &GC_zz,               NoDef_double,  NoMin_double,     NoMax_double      );
  
   ReadPara->Add( "GC_MASS",               &GC_mm,               NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "GC_RADIUS",               &GC_rr,               NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "GC_ANGLE",                &GC_theta,                0.0,         0.0,                   360.0      );
   

   ReadPara->Add( "FIX_CENTER",              &FixCenter,             Useless_bool,  Useless_bool,     Useless_bool      );
   ReadPara->Add( "SEARCH_RADIUS",           &SearchRadius,        NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "previous_center_x",           &previous_center_x,        NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "previous_center_y",           &previous_center_y,        NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "previous_center_z",           &previous_center_z,        NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Output_Profile_DT",           &output_profile_dt,        NoDef_double,  NoMin_double,     NoMax_double      );

   ReadPara->Read( FileName );
   delete ReadPara;
// (1-3) check the runtime parameters
   if ( OPT__INIT == INIT_BY_FUNCTION )
       Aux_Error( ERROR_INFO, "OPT__INIT=1 is not supported for this test problem !!\n" );
// (4) make a note
   if ( MPI_Rank == 0 )
      {
	Aux_Message( stdout, "=============================================================================\n" );
	Aux_Message( stdout, "  test problem ID = %d\n", TESTPROB_ID         );
	Aux_Message( stdout, "=============================================================================\n" );
      }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );



} // FUNCTION : SetParameter

//-------------------------------------------------------------------------------------------------------
// Function    :  AdjustGCPotential()
// Description :  Adjust the Potential
//
// Note        :  
// Parameter   :  add : the booldean value that determines the action (subtraction or add back)
//                GC_x,GC_y,GC_z : the current GC's position
//-------------------------------------------------------------------------------------------------------
void AdjustGCPotential(const bool add, const double GC_x, const double GC_y, const double GC_z)
{

   double sign = add ? 1.0 : -1.0;

   for (int lv=0; lv<NLEVEL; lv++)
   {
      const double dh = amr->dh[lv];
#  pragma omp for schedule( runtime )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         const double *EdgeL              = amr->patch[0][lv][PID]->EdgeL;
         real   (*PotPtr)[PS1][PS1] = amr->patch[ amr->PotSg[lv] ][lv][PID]->pot;
         const double x0                  = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
         const double y0                  = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
         const double z0                  = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;
         for (int k=0; k<PS1; k++)
         {
            const double z  = z0 + k*dh;
            const double dz = z - GC_z;
            for (int j=0; j<PS1; j++)
            {
               const double y  = y0 + j*dh;
               const double dy = y - GC_y;
               for (int i=0; i<PS1; i++)
               {
                  const double x  = x0 + i*dh;
                  const double dx = x - GC_x;
                  //const double r      = sqrt( SQR(dx) + SQR(dy) + SQR(dz) + SQR(dh) );
                  const double r      = sqrt( SQR(dx) + SQR(dy) + SQR(dz) );
                  const double GC_pot = -NEWTON_G*GC_mm/r;
         
                  PotPtr[k][j][i] += sign * (real)GC_pot;
               }
            }
          } 
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // for (int lv=0; lv<NLEVEL; lv++)
}



double interpolate(const double* x, const double* y, int size, double xi) {
    auto it = std::lower_bound(x, x + size, xi);
    if (it == x + size) return y[size - 1];  // xi is larger than any element in x
    if (it == x) return y[0];  // xi is smaller than any element in x

    int idx = std::distance(x, it);
    double x1 = x[idx - 1], x2 = x[idx];
    double y1 = y[idx - 1], y2 = y[idx];

    return y1 + (xi - x1) * (y2 - y1) / (x2 - x1);

}




double CalculateEnclosedMass(const double* radius, int radius_size,
                             const double* density, const double* weight,
                             double user_input_radius){

// Allocate arrays dynamically
//double* RADIUS = new double[radius_size];
//double* r_mid = new double[radius_size + 1];
//double* volume = new double[radius_size];
double* mass = new double[radius_size];
double* cumulative_mass = new double[radius_size];

double sum_mass = 0;
for (int i=0; i<radius_size;++i){
	mass[i] = density[i] * weight[i];
	sum_mass += mass[i];
	cumulative_mass[i] = sum_mass;
}




//double interpolated_mass_new = interpolate(r_mid, cumulative_mass, radius_size + 1, user_input_radius);
double interpolated_mass_new = interpolate(radius, cumulative_mass, radius_size, user_input_radius);

//if (MPI_Rank == 0 ) Aux_Message(stdout, "Interpolated mass at initial radius is %15.5e\n",interpolated_mass);
if (MPI_Rank == 0 ) Aux_Message(stdout, "New Interpolated mass at initial radius is %15.5e\n",interpolated_mass_new);

if ( MPI_Rank == 0 ){
      char Filename[MAX_STRING];
      sprintf( Filename, "Enclosed_Mass.txt" );
      FILE *File = fopen( Filename, "w" );
      fprintf( File, "#%19s  %21s  %21s \n", "Radius", "density", "enclosed_mass" );
      for (int b=0; b<radius_size; b++)
         fprintf( File, "%20.14e  %21.14e  %21.14e \n",
                  radius[b], density[b], cumulative_mass[b]  );
      fclose( File );
   
}

//delete[] RADIUS;
//delete[] r_mid;
//delete[] volume;
delete[] mass;
delete[] cumulative_mass;

return interpolated_mass_new;
}
//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_User_GC
// Description :  Update the GC position and Halo position in each step
//
// Note        :  1. Invoked by main() using the function pointer "Aux_Record_User_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__RECORD_USER"
//                3. This function will be called both during the program initialization and after each full update
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Aux_Record_User_GC()
{
   
// A. Find the GC's position
   double GC_x_local = 0.0, GC_y_local = 0.0, GC_z_local = 0.0;
   double GC_x, GC_y, GC_z;



   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++) {
      if ( amr->Par->Mass[p] < (real)0.0 )  continue;
      if ( amr->Par->Type[p] == PTYPE_GC)
      {
         GC_x_local = amr->Par->PosX[p];
         GC_y_local = amr->Par->PosY[p];
         GC_z_local = amr->Par->PosZ[p];

         char Filename[MAX_STRING];
         sprintf( Filename, "%s", "Record__Result_GC_position" );
         FILE *File = fopen( Filename, "a" );
         if ( Time[0]==0.0 )
         {
           fprintf(File, "%15s\t%15s\t%15s\t%15s\n", "#          Time", " GCX", " GCY", " GCZ");
         }
         fprintf(File, "%15.7f\t%15.7e\t%15.7e\t%15.7e\n", Time[0],GC_x_local,GC_y_local,GC_z_local );
         fclose ( File );
      }
   }

// A. Broadcast the GC's position to all rank
   MPI_Allreduce(&GC_x_local, &GC_x, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   MPI_Allreduce(&GC_y_local, &GC_y, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   MPI_Allreduce(&GC_z_local, &GC_z, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   
   MPI_Bcast(&GC_x,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&GC_y,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&GC_z,1,MPI_DOUBLE,0,MPI_COMM_WORLD);   
   MPI_Barrier(MPI_COMM_WORLD);

//   Aux_Message(stdout, "[%02d] GC's position is %15.6e %15.6e %15.6e\n",MPI_Rank,GC_x,GC_y,GC_z);



   // B. Set the center finding method base on Input__TestProb
   if ( FixCenter )
   {
      if ( MPI_Rank == 0 )
         {
         char Filename[MAX_STRING];
         sprintf( Filename, "%s", "Record__Result_Center" );
         FILE *File = fopen( Filename, "a" );
         if ( Time[0]==0.0 )
         {
            fprintf(File, "%15s\t%15s\t%15s\t%15s\n", "#          Time", " CenterX", " CenterY", " CenterZ");
         }
         fprintf(File, "%15.7f\t%15.7e\t%15.7e\t%15.7e\n",Time[0]  ,amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2]);
         fclose ( File );
      }
   }
   else // FixCenter is false
   {


   // B-1. subtract GC potential
      AdjustGCPotential(false,GC_x,GC_y,GC_z);     
      
       
   // B-2. Find the minimum potential position
      Extrema_t Extrema;
      Extrema.Field     = _POTE;
      Extrema.Radius    = SearchRadius*amr->dh[MAX_LEVEL];

      if ( Time[0]==0.0 )
      {
         Extrema.Center[0] = amr->BoxCenter[0];
         Extrema.Center[1] = amr->BoxCenter[1];
         Extrema.Center[2] = amr->BoxCenter[2];
         min_pot_last[0] = amr->BoxCenter[0];
         min_pot_last[1] = amr->BoxCenter[1];
         min_pot_last[2] = amr->BoxCenter[2];
      }
      else
      {
// Check if it is from RESTART

         if ( OPT__INIT == INIT_BY_RESTART )
         {
            Extrema.Center[0] = previous_center_x;
            Extrema.Center[1] = previous_center_y;
            Extrema.Center[2] = previous_center_z;
          
         }
         else
         {
         Extrema.Center[0] = min_pot_last[0];
         Extrema.Center[1] = min_pot_last[1];
         Extrema.Center[2] = min_pot_last[2];
         }      
      }

      Aux_FindExtrema( &Extrema, EXTREMA_MIN, 0, TOP_LEVEL, PATCH_LEAF );

      min_pot_last[0] = Extrema.Coord[0];
      min_pot_last[1] = Extrema.Coord[1];
      min_pot_last[2] = Extrema.Coord[2];

   // B-3. write them into a file

      if ( MPI_Rank == 0 )
      {
         char Filename[MAX_STRING];
         sprintf( Filename, "%s", "Record__Result_Center" );
         FILE *File = fopen( Filename, "a" );

         if (Time[0]==0.0)
         {
            fprintf(File, "%15s\t%15s\t%15s\t%15s\n", "#          Time", " CenterX", " CenterY", " CenterZ");
         }
         fprintf(File, "%15.7f\t%15.7e\t%15.7e\t%15.7e\n",
            Time[0]  ,Extrema.Coord[0], Extrema.Coord[1], Extrema.Coord[2]);
         fclose ( File );
      }

    // B-4 Add back the GC's potential
      AdjustGCPotential(true,GC_x,GC_y,GC_z);     


    // C. Output the profile if the step is 0
    //if ( Time[0] == 0.0 )
    if ( fabs(Time[0] - ideal_output_t) <= dTime_AllLv[0] )
    {
        // update the ideal output t
    //    Aux_Message(stdout,"[%02d] ideal_output_t is  %21.5f, dt is %21.5f \n",MPI_Rank, ideal_output_t, dTime_AllLv[0]);
        ideal_output_t += output_profile_dt;
         
        const double      Center[3]      = { Extrema.Coord[0],Extrema.Coord[1],Extrema.Coord[2] };
        //const double      Center[3]      = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[0]};
        const double      MaxRadius      = 0.6*amr->BoxSize[0]; // 0.6 to make sure every cell is counted
        const double      MinBinSize     = amr->dh[MAX_LEVEL];
        const bool        LogBin         = true;
        //const double      LogBinRatio    = 1.03;
        const double      LogBinRatio    = 1.0005;
        const bool        RemoveEmptyBin = true;
        const long        TVar[]         = { _DENS };
        const int         NProf          = 1;
        const int         MinLv          = 0;
        const int         MaxLv          = MAX_LEVEL;
        const PatchType_t PatchType      = PATCH_LEAF_PLUS_MAXNONLEAF;
        const double      PrepTime       = -1.0;

        Profile_t Prof_Dens;
        Profile_t *Prof[] = { &Prof_Dens };


        Aux_ComputeProfile( Prof, Center, MaxRadius, MinBinSize, LogBin, LogBinRatio, RemoveEmptyBin,
                            TVar, NProf, MinLv, MaxLv, PatchType, PrepTime );
        
        if ( MPI_Rank == 0 )
        {
           for (int p=0; p<NProf; p++)
           {
              char Filename[MAX_STRING];
              sprintf( Filename, "Density_Profile_%5f.txt",Time[0] );
              FILE *File = fopen( Filename, "w" );
              fprintf( File, "#%19s  %21s  %21s  %10s\n", "Radius", "Data", "Weight", "Cells" );
              for (int b=0; b<Prof[p]->NBin; b++)
                 fprintf( File, "%20.14e  %21.14e  %21.14e  %10ld\n",
                          Prof[p]->Radius[b], Prof[p]->Data[b], Prof[p]->Weight[b], Prof[p]->NCell[b] );
              fclose( File );
           }
        }
    }


   };

} // FUNCTION : Aux_Record_User_GC

#endif // #if ( MODEL == HYDRO )





//-------------------------------------------------------------------------------------------------------
// Function    :  Init_User_moveGC
// Description :  Template of user-defined initialization
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Init_User_Ptr",
//                   which must be set by a test problem initializer
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_User_moveGC()
{

if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

const long NNewPar = ( MPI_Rank == 0 ) ? 1 : 0;
const long NPar_AllRank = NNewPar;
real_par *NewParAttFlt[PAR_NATT_FLT_TOTAL];
long_par *NewParAttFlt[PAR_NATT_INT_TOTAL];

for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)   NewParAttFlt[v] = new real_par [NNewPar];
for (int v=0; v<PAR_NATT_INT_TOTAL; v++)   NewParAttInt[v] = new long_par [NNewPar];

// set particle attributes
real_par *Time_AllRank   =   NewParAttFlt[PAR_TIME];
real_par *Mass_AllRank   =   NewParAttFlt[PAR_MASS];
real_par *Pos_AllRank[3] = { NewParAttFlt[PAR_POSX], NewParAttFlt[PAR_POSY], NewParAttFlt[PAR_POSZ] };
real_par *Vel_AllRank[3] = { NewParAttFlt[PAR_VELX], NewParAttFlt[PAR_VELY], NewParAttFlt[PAR_VELZ] };
long_par *Type_AllRank   =   NewParAttInt[PAR_TYPE];
                


// find the center (base on maximum density)

Extrema_t Extrema;
Extrema.Field     = _DENS;
// Extrema.Radius = SearchRadius*amr->dh[MAX_LEVEL]; //the cell width of the finest level of resolution in the AMR grid multiply by the searchRadius setting
Extrema.Radius = __FLT_MAX__; // Entire Domain


Extrema.Center[0] = amr->BoxCenter[0];
Extrema.Center[1] = amr->BoxCenter[1];
Extrema.Center[2] = amr->BoxCenter[2];

Aux_FindExtrema( &Extrema, EXTREMA_MAX, 0, TOP_LEVEL, PATCH_LEAF );

Halo_Center_x_Initial = Extrema.Coord[0];
Halo_Center_y_Initial = Extrema.Coord[1];
Halo_Center_z_Initial = Extrema.Coord[2];

if ( MPI_Rank == 0 ){
Aux_Message( stdout, "-------------------------------------------\n");
Aux_Message( stdout, "Max Density x is at : %15.7e \n", Halo_Center_x_Initial);
Aux_Message( stdout, "Max Density y is at : %15.7e \n", Halo_Center_y_Initial);
Aux_Message( stdout, "Max Density z is at : %15.7e \n", Halo_Center_z_Initial);
Aux_Message( stdout, "-------------------------------------------\n");
}


// calculate the GC's velocity here
        double vc ;
        // update the ideal output t
        //    Aux_Message(stdout,"[%02d] ideal_output_t is  %21.5f, dt is %21.5f \n",MPI_Rank, ideal_output_t, dTime_AllLv[0]);
        // Set the velocity base on the calculation (only change the velocity)
         
        const double      Center[3]      = { Extrema.Coord[0],Extrema.Coord[1],Extrema.Coord[2] };
        const double      MaxRadius      = 0.6*amr->BoxSize[0]; // 0.6 to make sure every cell is counted
        const double      MinBinSize     = amr->dh[MAX_LEVEL];
        const bool        LogBin         = true;
        const double      LogBinRatio    = 1.0005;
        const bool        RemoveEmptyBin = true;
        const long        TVar[]         = { _DENS };
        const int         NProf          = 1;
        const int         MinLv          = 0;
        const int         MaxLv          = MAX_LEVEL;
        const PatchType_t PatchType      = PATCH_LEAF_PLUS_MAXNONLEAF;
        const double      PrepTime       = -1.0;

        Profile_t Prof_Dens;
        Profile_t *Prof[] = { &Prof_Dens };


        Aux_ComputeProfile( Prof, Center, MaxRadius, MinBinSize, LogBin, LogBinRatio, RemoveEmptyBin,
                            TVar, NProf, MinLv, MaxLv, PatchType, PrepTime );
        
        if ( OPT__INIT != INIT_BY_RESTART ) 
        {
        
            double* radius = Prof[0]->Radius;
            double* density = Prof[0]->Data;
            double* weight = Prof[0]->Weight;
            int radius_size = Prof[0]->NBin; // Number of bins

            
            double enclosed_mass = CalculateEnclosedMass(radius, radius_size, density,weight, GC_rr);
            
            
            vc = sqrt(NEWTON_G * enclosed_mass / GC_rr);
            
            double vc_relative = vc * enclosed_mass / (enclosed_mass + GC_mm);
        
        if ( MPI_Rank == 0 ) 
        {
           Aux_Message(stdout,"[%02d] velocity is %21.14e \n",MPI_Rank,vc);
        
           Aux_Message(stdout,"[%02d] enclosed mass is %21.14e \n",MPI_Rank,enclosed_mass);
           Aux_Message(stdout,"[%02d] radius is %21.14e \n",MPI_Rank,GC_rr);
           Aux_Message(stdout,"[%02d] Newton G is %21.14e \n",MPI_Rank,NEWTON_G);
           Aux_Message(stdout,"[%02d] GC Mass is %21.14e \n",MPI_Rank,GC_mm);
        }        
        
        
        }


    
   

   if ( MPI_Rank == 0 ){
      Aux_Message( stdout, "Setting GC's initial condition...\n");
      for (long p=0; p<NPar_AllRank; p++)
      {
         Time_AllRank[p] = Time[0];
         Mass_AllRank[p] = GC_mm;
         Type_AllRank[p] = PTYPE_GC;
         // Position
         Pos_AllRank[0][p] = Halo_Center_x_Initial + GC_rr*cos(GC_theta*M_PI/180);
         Pos_AllRank[1][p] = Halo_Center_y_Initial + GC_rr*sin(GC_theta*M_PI/180);
         Pos_AllRank[2][p] = Halo_Center_z_Initial;
         // Velocity
         Vel_AllRank[0][p] = -vc*sin(GC_theta*M_PI/180);
         Vel_AllRank[1][p] = vc*cos(GC_theta*M_PI/180);
         Vel_AllRank[2][p] = 0.0;
      } // for (long p=0; p<NPar_AllRank; p++)
      
      Aux_Message( stdout, "Setting GC's initial condition...done\n");

//    free memory
      delete [] ParData_AllRank;

   } // if ( MPI_Rank == 0)  


// add particles here
   Par_AddParticleAfterInit( NNewPar, NewParAttFlt, NewParAttInt );
   
// free memory
   for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)   delete [] NewParAttFlt[v];
   for (int v=0; v<PAR_NATT_INT_TOTAL; v++)   delete [] NewParAttInt[v];    


   if ( MPI_Rank == 0 )   Aux_Message( stdout, "%s ...done\n", __FUNCTION__ );

#  if ( defined PARTICLE  &&  defined LOAD_BALANCE )
   const double Par_Weight = amr->LB->Par_Weight;
#  else
   const double Par_Weight = 0.0;
#  endif
#  ifdef LOAD_BALANCE
   const UseLBFunc_t UseLB = USELB_YES;
#  else
   const UseLBFunc_t UseLB = USELB_NO;
#  endif

   for (int lv=0; lv<MAX_LEVEL; lv++)
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Refining level %d ...\n", lv );

      Flag_Real( lv, UseLB );

      Refine( lv, UseLB );

#     ifdef LOAD_BALANCE
      LB_Init_LoadBalance( true, true, Par_Weight, true, false, lv+1 );
#     endif

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Refining level %d ... done\n", lv );
   } // for (int lv=OPT__UM_IC_LEVEL; lv<MAX_LEVEL; lv++)


} // FUNCTION : Init_User_moveGC





//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_DynamicalFriction
// Description :  Test problem initializer for dynamical friction problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_Dynamical_Friction()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();

#  if ( MODEL == ELBDM && defined GRAVITY )
// set the problem-specific runtime parameters
   SetParameter();
   Aux_Record_User_Ptr     = Aux_Record_User_GC;
   Init_User_Ptr	   = Init_User_moveGC;
#  endif
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_GC
