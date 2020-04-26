/**
 * @file    PIAAACMCsimul.c
 * @brief   PIAA-type coronagraph design
 * 
 * Can design both APLCMC and PIAACMC coronagraphs
 *  
 *
 * 
 */



/* ================================================================== */
/* ================================================================== */
/*            MODULE INFO                                             */
/* ================================================================== */
/* ================================================================== */

// module default short name
// all CLI calls to this module functions will be <shortname>.<funcname>
// if set to "", then calls use <funcname>
#define MODULE_SHORTNAME_DEFAULT "piaacmcsim"

// Module short description
#define MODULE_DESCRIPTION       "PIAA complex mask coronagraph simulation"





/* ================================================================== */
/* ================================================================== */
/*            DEPENDANCIES                                            */
/* ================================================================== */
/* ================================================================== */


#define _GNU_SOURCE


// System include

#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>
#include <time.h>

// External libraries

#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <fitsio.h>


// milk includes
//   core modules
#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
//   other modules
#include "info/info.h"
#include "fft/fft.h"
#include "image_gen/image_gen.h"
#include "WFpropagate/WFpropagate.h"
#include "statistic/statistic.h"
#include "linopt_imtools/linopt_imtools.h"
#include "OpticsMaterials/OpticsMaterials.h"
#include "image_filter/image_filter.h"
#include "image_basic/image_basic.h"
#include "coronagraphs/coronagraphs.h"
#include "PIAACMCsimul/PIAACMCsimul.h"
#include "OptSystProp/OptSystProp.h"



# ifdef HAVE_LIBGOMP
#include <omp.h>
#define OMP_NELEMENT_LIMIT 1000000
#endif




/* ================================================================== */
/* ================================================================== */
/*            GLOBAL VARIABLES                                        */
/* ================================================================== */
/* ================================================================== */


static int INITSTATUS_PIAACMCsimul = 0;



PIAACMCsimul_varType piaacmcsimul_var; 

/// optical system description
OPTSYST *optsyst;


OPTPIAACMCDESIGN *piaacmc;





/* ================================================================== */
/* ================================================================== */
/*            INITIALIZE LIBRARY                                      */
/* ================================================================== */
/* ================================================================== */

// Module initialization macro in CLIcore.h
// macro argument defines module name for bindings
//
INIT_MODULE_LIB(milk_module_example)







/* ================================================================== */
/* ================================================================== */
/*            COMMAND LINE INTERFACE (CLI) FUNCTIONS                  */
/* ================================================================== */
/* ================================================================== */





// command line interface (CLI) commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string, not existing image
// 4: existing image
// 5: string 
//



// local function(s)
static double f_evalmask (const gsl_vector *v, void *params);

// for testing
static char flogcomment[200];
//#define PIAASIMUL_LOGFUNC0 // top level
//#define PIAASIMUL_LOGFUNC1 // lower level



/* =============================================================================================== */
/* =============================================================================================== */
/** @name  Command line interface (CLI)
 *  CLI commands
 */
///@{
/* =============================================================================================== */
/* =============================================================================================== */


/* =============================================================================================== */
/*  2. Focal plane mask construction                                                               */
/* =============================================================================================== */

errno_t PIAACMCsimul_rings2sectors__cli()
{
    if(0
            + CLI_checkarg(1, CLIARG_IMG)
            + CLI_checkarg(2, CLIARG_STR_NOT_IMG)
            + CLI_checkarg(3, CLIARG_STR_NOT_IMG)
            == 0)
    {
        PIAACMCsimul_rings2sectors(
            data.cmdargtoken[1].val.string,
            data.cmdargtoken[2].val.string,
            data.cmdargtoken[3].val.string);

        return CLICMD_SUCCESS;
    }
    else
    {
        return CLICMD_INVALID_ARG;
    }
}


/* =============================================================================================== */
/*  4. Lyot Stop(s)                                                              */
/* =============================================================================================== */

errno_t PIAACMCsimul_geomProp__cli()
{
    if(0
            + CLI_checkarg(1, CLIARG_IMG)
            + CLI_checkarg(2, CLIARG_IMG)
            + CLI_checkarg(3, CLIARG_STR_NOT_IMG)
            + CLI_checkarg(4, CLIARG_STR_NOT_IMG)
            + CLI_checkarg(5, CLIARG_FLOAT)
            + CLI_checkarg(6, CLIARG_FLOAT)
            + CLI_checkarg(7, CLIARG_FLOAT)
            + CLI_checkarg(8, CLIARG_FLOAT)
            + CLI_checkarg(9, CLIARG_FLOAT)
            + CLI_checkarg(10, CLIARG_FLOAT)
            == 0)
    {
        PIAACMCsimul_geomProp(
            data.cmdargtoken[1].val.string,
            data.cmdargtoken[2].val.string,
            data.cmdargtoken[3].val.string,
            data.cmdargtoken[4].val.string,
            data.cmdargtoken[5].val.numf,
            data.cmdargtoken[6].val.numf,
            data.cmdargtoken[7].val.numf,
            data.cmdargtoken[8].val.numf,
            data.cmdargtoken[9].val.numf,
            data.cmdargtoken[10].val.numf);

        return CLICMD_SUCCESS;
    }
    else
    {
        return CLICMD_INVALID_ARG;
    }

}


/* =============================================================================================== */
/*  5. Focal plane mask optimization                                                               */
/* =============================================================================================== */


errno_t PIAACMC_FPMresp_rmzones__cli()
{
    if(0
            + CLI_checkarg(1, CLIARG_IMG)
            + CLI_checkarg(2, CLIARG_STR_NOT_IMG)
            + CLI_checkarg(3, CLIARG_LONG)
            == 0)
    {
        PIAACMC_FPMresp_rmzones(
            data.cmdargtoken[1].val.string,
            data.cmdargtoken[2].val.string,
            data.cmdargtoken[3].val.numl);

        return CLICMD_SUCCESS;
    }
    else
    {
        return CLICMD_INVALID_ARG;
    }
}


errno_t PIAACMC_FPMresp_resample__cli()
{
    if(0
            + CLI_checkarg(1, CLIARG_IMG)
            + CLI_checkarg(2, CLIARG_STR_NOT_IMG)
            + CLI_checkarg(3, CLIARG_LONG)
            + CLI_checkarg(4, CLIARG_LONG)
            == 0)
    {
        PIAACMC_FPMresp_resample(
            data.cmdargtoken[1].val.string,
            data.cmdargtoken[2].val.string,
            data.cmdargtoken[3].val.numl,
            data.cmdargtoken[4].val.numl);
            
        return CLICMD_SUCCESS;
    }
    else
    {
        return CLICMD_INVALID_ARG;
    }
}


/* =============================================================================================== */
/*  6. Focal plane processing                                                                      */
/* =============================================================================================== */

errno_t PIAACMC_FPM_process__cli()
{
    if(0
            + CLI_checkarg(1, CLIARG_IMG)
            + CLI_checkarg(2, CLIARG_STR)
            + CLI_checkarg(3, CLIARG_LONG)
            + CLI_checkarg(4, CLIARG_STR_NOT_IMG)
            == 0)
    {
        PIAACMC_FPM_process(
            data.cmdargtoken[1].val.string,
            data.cmdargtoken[2].val.string,
            data.cmdargtoken[3].val.numl,
            data.cmdargtoken[4].val.string);

        return CLICMD_SUCCESS;
    }
    else
    {
        return CLICMD_INVALID_ARG;
    }
}


/* =============================================================================================== */
/*  7. High level routines                                                                         */
/* =============================================================================================== */

errno_t PIAACMCsimul_run__cli()
{
    if(0
            + CLI_checkarg(1, CLIARG_STR_NOT_IMG)
            + CLI_checkarg(2, CLIARG_LONG)
            == 0)
    {
        PIAACMCsimul_run(
            data.cmdargtoken[1].val.string,
            data.cmdargtoken[2].val.numl);

        return CLICMD_SUCCESS;
    }
    else
    {
        return CLICMD_INVALID_ARG;
    }
}

///@}





 /** @name MODULE INITIALIZATION
  * Registers CLI commands
 */
///@{
 
static errno_t init_module_CLI()
{

    RegisterCLIcommand(
        "piaacmcsimring2sect",
        __FILE__,
        PIAACMCsimul_rings2sectors__cli,
        "turn ring fpm design into sectors",
        "<input ring fpm> <zone-ring table> <output sector fpm>",
        "piaacmcsimring2sect",
        "long PIAACMCsimul_rings2sectors(const char *IDin_name, const char *sectfname, const char *IDout_name)");

	
    RegisterCLIcommand(
        "piaacmcsimrun",
        __FILE__,
        PIAACMCsimul_run__cli,
        "Simulate PIAACMC",
        "<configuration index [string]> <mode[int]>",
        "piaacmcsimrun",
        "int PIAACMCsimul_run(const char *confindex, long mode)");

  	
    RegisterCLIcommand(
        "piaacmsimfpmresprm",
        __FILE__,
        PIAACMC_FPMresp_rmzones__cli,
        "remove zones in FPM resp matrix",
        "<input FPMresp> <output FPMresp> <NBzone removed>",
        "piaacmsimfpmresprm FPMresp FPMrespout 125",
        "long PIAACMC_FPMresp_rmzones(const char *FPMresp_in_name, const char *FPMresp_out_name, long NBzones)");


    RegisterCLIcommand(
        "piaacmsimfpmresprs",
        __FILE__,
        PIAACMC_FPMresp_resample__cli,
        "resample FPM resp matrix",
        "<input FPMresp> <output FPMresp> <NBlambda> <EvalPts step>",
        "piaacmsimfpmresprs FPMresp FPMrespout 10 2",
        "long PIAACMC_FPMresp_resample(const char *FPMresp_in_name, const char *FPMresp_out_name, long NBlambda, long PTstep)");


    RegisterCLIcommand(
        "piaacmcfpmprocess",
        __FILE__,
        PIAACMC_FPM_process__cli,
        "Quantize FPM",
        "<input FPM sags> <sectors ASCII file> <number of exposures> <output FPM sags>",
        "piaacmcfpmprocess",
        "long PIAACMC_FPM_process(const char *FPMsag_name, const char *zonescoord_name, long NBexp, const char *outname)");

    RegisterCLIcommand(
        "piaacmcgeomprop",
        __FILE__,
        PIAACMCsimul_geomProp__cli,
        "Geometric propagation from surface",
        "<input map> <input sag> <output map> <output map counter> <delta refractive index> <pixel scale> <propagation dist> <kernel radius> <kernel step> <clear aperture>",
        "piaacmcgeomprop pupin piaa0z pupout cntout 2.0 0.00011 2.302606 3.0 0.5 200.0",
        "long PIAACMCsimul_geomProp(const char *IDin_name, const char *IDsag_name, const char *IDout_name, const char *IDoutcnt_name, float drindex, float pscale, float zprop, float krad, float kstep, float rlim)");
    
    


	piaacmcsimul_var.optsystinit = 0;
	
	// this makes 20% bandwidth, from 0.55/1.1 to 0.55*1.1
	piaacmcsimul_var.LAMBDASTART = 0.5e-6;
	piaacmcsimul_var.LAMBDAEND = 0.605e-6;

	piaacmcsimul_var.FORCE_CREATE_Cmodes = 0;
	piaacmcsimul_var.CREATE_Cmodes = 0;
	piaacmcsimul_var.FORCE_CREATE_Fmodes = 0;
	piaacmcsimul_var.CREATE_Fmodes = 0;

	piaacmcsimul_var.FORCE_CREATE_fpmzmap = 0;
	piaacmcsimul_var.CREATE_fpmzmap = 0;
	piaacmcsimul_var.FORCE_CREATE_fpmzt = 0;
	piaacmcsimul_var.CREATE_fpmzt = 0;

	piaacmcsimul_var.FORCE_CREATE_fpmza = 0;
	piaacmcsimul_var.CREATE_fpmza;

	piaacmcsimul_var.FORCE_MAKE_PIAA0shape = 0;
	piaacmcsimul_var.MAKE_PIAA0shape = 0;
	piaacmcsimul_var.FORCE_MAKE_PIAA1shape = 0;
	piaacmcsimul_var.MAKE_PIAA1shape = 0;
	
	piaacmcsimul_var.focmMode = -1; // if != -1, compute only impulse response to corresponding zone
	piaacmcsimul_var.PIAACMC_FPMsectors = 0;

	piaacmcsimul_var.FPMSCALEFACTOR = 0.9;
	piaacmcsimul_var.PIAACMC_MASKRADLD = 0.0;
	
	piaacmcsimul_var.LOOPCNT = 0;
	
	piaacmcsimul_var.CnormFactor = 1.0;

	piaacmcsimul_var.computePSF_FAST_FPMresp = 0;
	piaacmcsimul_var.computePSF_ResolvedTarget = 0; // source size = 1e-{0.1*computePSF_ResolvedTarget}
	piaacmcsimul_var.computePSF_ResolvedTarget_mode = 0; // 0: source is simulated as 3 points, 1: source is simulated as 6 points
	piaacmcsimul_var.PIAACMC_FPM_FASTDERIVATIVES = 0;

	piaacmcsimul_var.SCORINGTOTAL = 1.0;
	piaacmcsimul_var.MODampl = 1.0e-6;
	piaacmcsimul_var.SCORINGMASKTYPE = 0;
	piaacmcsimul_var.PIAACMC_save = 1;
	//piaacmcsimul_var.PIAACMC_MASKregcoeff = 1.0;
	piaacmcsimul_var.PIAACMC_fpmtype = 0;
	
	piaacmcsimul_var.WRITE_OK = 1;
	
	piaacmcsimul_var.LINOPT = 0;
	
    // add atexit functions here
    atexit(PIAACMCsimul_free);

    return RETURN_SUCCESS;

}


///@}


// first argument should be "PIAACMCsimul.fcall.log"
// second argument should be __FUNCTION__
// PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, "");
static void PIAACMCsimul_logFunctionCall(
    char *LogFileName,
    const char *FunctionName,
    long line,
    char *comments
)
{
    FILE *fp;
    time_t tnow;
    struct tm *uttime;
    struct timespec timenow;

    char string[21];

    tnow = time(NULL);
    uttime = gmtime(&tnow);
    clock_gettime(CLOCK_REALTIME, &timenow);

    // add custom parameter
    if(piaacmc == NULL)
    {
        sprintf(string, "NULL");
    }
    else
    {
        sprintf(string, "%20ld", piaacmc[0].focmNBzone);
    }

    fp = fopen(LogFileName, "a");
    fprintf(fp, "%02d:%02d:%02ld.%09ld  %10d  %40s %6ld   %20s %s\n",
            uttime->tm_hour, uttime->tm_min, timenow.tv_sec % 60, timenow.tv_nsec, getpid(),
            FunctionName, line, string, comments);
    fclose(fp);
}





















































































































