#include "GAMER.h"
#include "TestProb.h"



static void OutputError();


// problem-specific global variables
// =======================================================================================
static double Acoustic_RhoAmp;      // amplitude of the density perturbation (assuming background density = 1.0)
static double Acoustic_Cs;          // sound speed
static double Acoustic_v0;          // background velocity
static int    Acoustic_Dir;         // wave direction: (0/1/2/3) --> (x/y/z/diagonal)
static double Acoustic_Sign;        // (+1/-1) --> (right/left-moving wave)
static double Acoustic_Phase0;      // initial phase shift

static double Acoustic_WaveLength;  // wavelength
// =======================================================================================




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
#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

#  ifndef FLOAT8
   Aux_Error( ERROR_INFO, "FLOAT8 must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

#  if ( EOS != EOS_GAMMA )
   Aux_Error( ERROR_INFO, "EOS != EOS_GAMMA !!\n" );
#  endif

   if ( Acoustic_Dir == 3  &&  ( amr->BoxSize[0] != amr->BoxSize[1] || amr->BoxSize[0] != amr->BoxSize[2] )  )
      Aux_Error( ERROR_INFO, "simulation domain must be cubic for Acoustic_Dir = %d !!\n", Acoustic_Dir );

   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] != BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "please set \"OPT__BC_FLU_* = 1\" (i.e., periodic BC) !!\n" );


// warnings
   if ( MPI_Rank == 0 )
   {
      if ( !OPT__OUTPUT_USER )   Aux_Message( stdout, "WARNING : OPT__OUTPUT_USER is off !!\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO )
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

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Acoustic_RhoAmp",   &Acoustic_RhoAmp,       -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Acoustic_Cs",       &Acoustic_Cs,           -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Acoustic_v0",       &Acoustic_v0,            0.0,          NoMin_double,     NoMax_double      );
   ReadPara->Add( "Acoustic_Dir",      &Acoustic_Dir,           3,            0,                3                 );
   ReadPara->Add( "Acoustic_Sign",     &Acoustic_Sign,          1.0,          NoMin_double,     NoMax_double      );
   ReadPara->Add( "Acoustic_Phase0",   &Acoustic_Phase0,        0.0,          NoMin_double,     NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;


// (2) set the problem-specific derived parameters
   Acoustic_WaveLength = ( Acoustic_Dir == 3 ) ? amr->BoxSize[0]/sqrt(3.0) : amr->BoxSize[Acoustic_Dir];

// force Acoustic_Sign to be +1.0/-1.0
   if ( Acoustic_Sign >= 0.0 )   Acoustic_Sign = +1.0;
   else                          Acoustic_Sign = -1.0;


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const double End_T_Default    = Acoustic_WaveLength / Acoustic_Cs;
   const long   End_Step_Default = __INT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID     = %d\n",      TESTPROB_ID );
      Aux_Message( stdout, "  density amplitude   = % 14.7e\n", Acoustic_RhoAmp );
      Aux_Message( stdout, "  sound speed         = % 14.7e\n", Acoustic_Cs );
      Aux_Message( stdout, "  background velocity = % 14.7e\n", Acoustic_v0 );
      Aux_Message( stdout, "  direction           = %d\n",      Acoustic_Dir );
      Aux_Message( stdout, "  sign (R/L)          = % 14.7e\n", Acoustic_Sign );
      Aux_Message( stdout, "  initial phase shift = % 14.7e\n", Acoustic_Phase0 );
      Aux_Message( stdout, "  wavelength          = % 14.7e\n", Acoustic_WaveLength );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

// assuming EOS_GAMMA
   const double _Gamma_m1 = 1.0/(GAMMA-1.0);

   double r, v1, P0, P1, Phase, WaveK, WaveW;
   double Dens, Mom, MomX, MomY, MomZ, Pres, Eint, Etot;

   switch ( Acoustic_Dir ) {
      case 0:  r = x;                        break;
      case 1:  r = y;                        break;
      case 2:  r = z;                        break;
      case 3:  r = ( x + y + z )/sqrt(3.0);  break;
   }
   r -= Acoustic_v0*Time;

   v1    = Acoustic_Sign*Acoustic_Cs*Acoustic_RhoAmp;
   P0    = SQR(Acoustic_Cs)/GAMMA;
   P1    = SQR(Acoustic_Cs)*Acoustic_RhoAmp;

   WaveK = 2.0*M_PI/Acoustic_WaveLength;
   WaveW = 2.0*M_PI/(Acoustic_WaveLength/Acoustic_Cs);
   Phase = WaveK*r - Acoustic_Sign*WaveW*Time + Acoustic_Phase0;

   Dens  = 1.0 + Acoustic_RhoAmp*cos(Phase);
   Mom   = Dens*( v1*cos(Phase) + Acoustic_v0 );

   switch ( Acoustic_Dir ) {
      case 0:  MomX = Mom;            MomY = 0.0;   MomZ = 0.0;   break;
      case 1:  MomX = 0.0;            MomY = Mom;   MomZ = 0.0;   break;
      case 2:  MomX = 0.0;            MomY = 0.0;   MomZ = Mom;   break;
      case 3:  MomX = Mom/sqrt(3.0);  MomY = MomX;  MomZ = MomX;  break;
   }

   Pres  = P0 + P1*cos(Phase);
   Eint  = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                     EoS_AuxArray_Int, h_EoS_Table );   // assuming EoS requires no passive scalars
   Etot  = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here

// set the output array
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  OutputError
// Description :  Output the L1 error
//
// Note        :  1. Invoke Output_L1Error()
//                2. Use SetGridIC() to provide the analytical solution at any given time
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void OutputError()
{

   const char Prefix[100]     = "AcousticWave";
   const OptOutputPart_t Part = OUTPUT_X + Acoustic_Dir;

   Output_L1Error( SetGridIC, NULL, Prefix, Part, OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z );

} // FUNCTION : OutputError
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_AcousticWave
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_AcousticWave()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr = SetGridIC;
   Output_User_Ptr        = OutputError;
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_AcousticWave
