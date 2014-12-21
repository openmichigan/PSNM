/*
 * clInfo.c --
 *
 *      Program to enumerate and dump all of the OpenCL information for a
 *      machine (or at least for a specific run-time).
 *
 * (Jeremy Sugerman, 13 August 2009)
 */

#include <stdio.h>
#include <malloc.h>
#include <getopt.h>
#include <stdlib.h>

#include "CL/cl.h"

#define Warning(...)    fprintf(stderr, __VA_ARGS__)

typedef struct Opts {
   /*
    * Boolean options below here, all default to false (zero).
    */

   int verify;
   int timing;
} Opts;


/*
 * Usage --
 *
 *      Prints the usage message and exits.
 *
 * Results:
 *      void, but calls exit(1)...
 */

static void
Usage(const char *progName) {
   Warning("Usage: %s [options]\n", progName);
   Warning("Options:\n");
   Warning("  -h, --help                This message\n");

   exit(1);
}


/*
 * ParseOpts --
 *
 *      Converts the commandline parameters into their internal
 *      representation.
 *
 * Results:
 *      void, opts is initialized.
 */

static void
ParseOpts(Opts *opts, int argc, char *argv[]) {
   int opt;

   static struct option longOptions[] = {
      {"help",         0, 0, 'h'},
   };


   while ((opt = getopt_long(argc, argv, "h",
                             longOptions, NULL)) != EOF) {
      switch(opt) {
      case 'h':
      default:
         Usage(argv[0]);
         break;
      }
   }

   return;
}


/*
 * CLErrString --
 *
 *      Utility function that converts an OpenCL status into a human
 *      readable string.
 *
 * Results:
 *      const char * pointer to a static string.
 */

static const char *
CLErrString(cl_int status) {
   static struct { cl_int code; const char *msg; } error_table[] = {
      { CL_SUCCESS, "success" },
      { CL_DEVICE_NOT_FOUND, "device not found", },
      { CL_DEVICE_NOT_AVAILABLE, "device not available", },
      { CL_COMPILER_NOT_AVAILABLE, "compiler not available", },
      { CL_MEM_OBJECT_ALLOCATION_FAILURE, "mem object allocation failure", },
      { CL_OUT_OF_RESOURCES, "out of resources", },
      { CL_OUT_OF_HOST_MEMORY, "out of host memory", },
      { CL_PROFILING_INFO_NOT_AVAILABLE, "profiling not available", },
      { CL_MEM_COPY_OVERLAP, "memcopy overlaps", },
      { CL_IMAGE_FORMAT_MISMATCH, "image format mismatch", },
      { CL_IMAGE_FORMAT_NOT_SUPPORTED, "image format not supported", },
      { CL_BUILD_PROGRAM_FAILURE, "build program failed", },
      { CL_MAP_FAILURE, "map failed", },
      { CL_INVALID_VALUE, "invalid value", },
      { CL_INVALID_DEVICE_TYPE, "invalid device type", },
      { 0, NULL },
   };
   static char unknown[25];
   int ii;

   for (ii = 0; error_table[ii].msg != NULL; ii++) {
      if (error_table[ii].code == status) {
         return error_table[ii].msg;
      }
   }

   snprintf(unknown, sizeof unknown, "unknown error %d", status);
   return unknown;
}


/*
 * PrintDevice --
 *
 *      Dumps everything about the given device ID.
 *
 * Results:
 *      void.
 */

static void
PrintDevice(cl_device_id device) {
#define LONG_PROPS \
  defn(VENDOR_ID), \
  defn(MAX_COMPUTE_UNITS), \
  defn(MAX_WORK_ITEM_DIMENSIONS), \
  defn(MAX_WORK_GROUP_SIZE), \
  defn(PREFERRED_VECTOR_WIDTH_CHAR), \
  defn(PREFERRED_VECTOR_WIDTH_SHORT), \
  defn(PREFERRED_VECTOR_WIDTH_INT), \
  defn(PREFERRED_VECTOR_WIDTH_LONG), \
  defn(PREFERRED_VECTOR_WIDTH_FLOAT), \
  defn(PREFERRED_VECTOR_WIDTH_DOUBLE), \
  defn(MAX_CLOCK_FREQUENCY), \
  defn(ADDRESS_BITS), \
  defn(MAX_MEM_ALLOC_SIZE), \
  defn(IMAGE_SUPPORT), \
  defn(MAX_READ_IMAGE_ARGS), \
  defn(MAX_WRITE_IMAGE_ARGS), \
  defn(IMAGE2D_MAX_WIDTH), \
  defn(IMAGE2D_MAX_HEIGHT), \
  defn(IMAGE3D_MAX_WIDTH), \
  defn(IMAGE3D_MAX_HEIGHT), \
  defn(IMAGE3D_MAX_DEPTH), \
  defn(MAX_SAMPLERS), \
  defn(MAX_PARAMETER_SIZE), \
  defn(MEM_BASE_ADDR_ALIGN), \
  defn(MIN_DATA_TYPE_ALIGN_SIZE), \
  defn(GLOBAL_MEM_CACHELINE_SIZE), \
  defn(GLOBAL_MEM_CACHE_SIZE), \
  defn(GLOBAL_MEM_SIZE), \
  defn(MAX_CONSTANT_BUFFER_SIZE), \
  defn(MAX_CONSTANT_ARGS), \
  defn(LOCAL_MEM_SIZE), \
  defn(ERROR_CORRECTION_SUPPORT), \
  defn(PROFILING_TIMER_RESOLUTION), \
  defn(ENDIAN_LITTLE), \
  defn(AVAILABLE), \
  defn(COMPILER_AVAILABLE),

#define STR_PROPS \
  defn(NAME), \
  defn(VENDOR), \
  defn(PROFILE), \
  defn(VERSION), \
  defn(EXTENSIONS),

#define HEX_PROPS \
   defn(SINGLE_FP_CONFIG), \
   defn(QUEUE_PROPERTIES),


/* XXX For completeness, it'd be nice to dump this one, too. */
#define WEIRD_PROPS \
   CL_DEVICE_MAX_WORK_ITEM_SIZES,

   static struct { cl_device_info param; const char *name; } longProps[] = {
#define defn(X) { CL_DEVICE_##X, #X }
      LONG_PROPS
#undef defn
      { 0, NULL },
   };
   static struct { cl_device_info param; const char *name; } hexProps[] = {
#define defn(X) { CL_DEVICE_##X, #X }
      HEX_PROPS
#undef defn
      { 0, NULL },
   };
   static struct { cl_device_info param; const char *name; } strProps[] = {
#define defn(X) { CL_DEVICE_##X, #X }
      STR_PROPS
#undef defn
      { CL_DRIVER_VERSION, "DRIVER_VERSION" },
      { 0, NULL },
   };
   cl_int status;
   size_t size;
   char buf[65536];
   long long val; /* Avoids unpleasant surprises for some params */
   int ii;


   for (ii = 0; strProps[ii].name != NULL; ii++) {
      status = clGetDeviceInfo(device, strProps[ii].param, sizeof buf, buf, &size);
      if (status != CL_SUCCESS) {
         Warning("\tdevice[%p]: Unable to get %s: %s!\n",
                 device, strProps[ii].name, CLErrString(status));
         continue;
      }
      if (size > sizeof buf) {
         Warning("\tdevice[%p]: Large %s (%d bytes)!  Truncating to %d!\n",
                 device, strProps[ii].name, size, sizeof buf);
      }
      printf("\tdevice[%p]: %s: %s\n",
             device, strProps[ii].name, buf);
   }
   printf("\n");

   status = clGetDeviceInfo(device, CL_DEVICE_TYPE, sizeof val, &val, NULL);
   if (status == CL_SUCCESS) {
      printf("\tdevice[%p]: Type: ", device);
      if (val & CL_DEVICE_TYPE_DEFAULT) {
         val &= ~CL_DEVICE_TYPE_DEFAULT;
         printf("Default ");
      }
      if (val & CL_DEVICE_TYPE_CPU) {
         val &= ~CL_DEVICE_TYPE_CPU;
         printf("CPU ");
      }
      if (val & CL_DEVICE_TYPE_GPU) {
         val &= ~CL_DEVICE_TYPE_GPU;
         printf("GPU ");
      }
      if (val & CL_DEVICE_TYPE_ACCELERATOR) {
         val &= ~CL_DEVICE_TYPE_ACCELERATOR;
         printf("Accelerator ");
      }
      if (val != 0) {
         printf("Unknown (0x%llx) ", val);
      }
      printf("\n");
   } else {
      Warning("\tdevice[%p]: Unable to get TYPE: %s!\n",
              device, CLErrString(status));
   }

   status = clGetDeviceInfo(device, CL_DEVICE_EXECUTION_CAPABILITIES,
                            sizeof val, &val, NULL);
   if (status == CL_SUCCESS) {
      printf("\tdevice[%p]: EXECUTION_CAPABILITIES: ", device);
      if (val & CL_EXEC_KERNEL) {
         val &= ~CL_EXEC_KERNEL;
         printf("Kernel ");
      }
      if (val & CL_EXEC_NATIVE_KERNEL) {
         val &= ~CL_EXEC_NATIVE_KERNEL;
         printf("Native ");
      }
      if (val) {
         printf("Unknown (0x%llx) ", val);
      }
      printf("\n");
   } else {
      Warning("\tdevice[%p]: Unable to get EXECUTION_CAPABILITIES: %s!\n",
              device, CLErrString(status));
   }

   status = clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_CACHE_TYPE,
                            sizeof val, &val, NULL);
   if (status == CL_SUCCESS) {
      static const char *cacheTypes[] = { "None", "Read-Only", "Read-Write" };
      static int numTypes = sizeof cacheTypes / sizeof cacheTypes[0];

      printf("\tdevice[%p]: GLOBAL_MEM_CACHE_TYPE: %s (%lld)\n",
             device, val < numTypes ? cacheTypes[val] : "???", val);
   } else {
      Warning("\tdevice[%p]: Unable to get GLOBAL_MEM_CACHE_TYPE: %s!\n",
              device, CLErrString(status));
   }
   status = clGetDeviceInfo(device,
                            CL_DEVICE_LOCAL_MEM_TYPE, sizeof val, &val, NULL);
   if (status == CL_SUCCESS) {
      static const char *lmemTypes[] = { "???", "Local", "Global" };
      static int numTypes = sizeof lmemTypes / sizeof lmemTypes[0];

      printf("\tdevice[%p]: CL_DEVICE_LOCAL_MEM_TYPE: %s (%lld)\n",
             device, val < numTypes ? lmemTypes[val] : "???", val);
   } else {
      Warning("\tdevice[%p]: Unable to get CL_DEVICE_LOCAL_MEM_TYPE: %s!\n",
              device, CLErrString(status));
   }

   for (ii = 0; hexProps[ii].name != NULL; ii++) {
      status = clGetDeviceInfo(device, hexProps[ii].param, sizeof val, &val, &size);
      if (status != CL_SUCCESS) {
         Warning("\tdevice[%p]: Unable to get %s: %s!\n",
                 device, hexProps[ii].name, CLErrString(status));
         continue;
      }
      if (size > sizeof val) {
         Warning("\tdevice[%p]: Large %s (%d bytes)!  Truncating to %d!\n",
                 device, hexProps[ii].name, size, sizeof val);
      }
      printf("\tdevice[%p]: %s: 0x%llx\n",
             device, hexProps[ii].name, val);
   }
   printf("\n");

   for (ii = 0; longProps[ii].name != NULL; ii++) {
      status = clGetDeviceInfo(device, longProps[ii].param, sizeof val, &val, &size);
      if (status != CL_SUCCESS) {
         Warning("\tdevice[%p]: Unable to get %s: %s!\n",
                 device, longProps[ii].name, CLErrString(status));
         continue;
      }
      if (size > sizeof val) {
         Warning("\tdevice[%p]: Large %s (%d bytes)!  Truncating to %d!\n",
                 device, longProps[ii].name, size, sizeof val);
      }
      printf("\tdevice[%p]: %s: %lld\n",
             device, longProps[ii].name, val);
   }
}


/*
 * PrintPlatform --
 *
 *      Dumps everything about the given platform ID.
 *
 * Results:
 *      void.
 */

static void
PrintPlatform(cl_platform_id platform) {
   static struct { cl_platform_info param; const char *name; } props[] = {
      { CL_PLATFORM_PROFILE, "profile" },
      { CL_PLATFORM_VERSION, "version" },
      { CL_PLATFORM_NAME, "name" },
      { CL_PLATFORM_VENDOR, "vendor" },
      { CL_PLATFORM_EXTENSIONS, "extensions" },
      { 0, NULL },
   };
   cl_device_id *deviceList;
   cl_uint numDevices;
   cl_int status;
   char buf[65536];
   size_t size;
   int ii;

   for (ii = 0; props[ii].name != NULL; ii++) {
      status = clGetPlatformInfo(platform, props[ii].param, sizeof buf, buf, &size);
      if (status != CL_SUCCESS) {
         Warning("platform[%p]: Unable to get %s: %s\n",
                 platform, props[ii].name, CLErrString(status));
         continue;
      }
      if (size > sizeof buf) {
         Warning("platform[%p]: Huge %s (%d bytes)!  Truncating to %d\n",
                 platform, props[ii].name, size, sizeof buf);
      }
      printf("platform[%p]: %s: %s\n", platform, props[ii].name, buf);
   }

   if ((status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL,
                                0, NULL, &numDevices)) != CL_SUCCESS) {
      Warning("platform[%p]: Unable to query the number of devices: %s\n",
              platform, CLErrString(status));
      return;
   }
   printf("platform[%p]: Found %d device(s).\n", platform, numDevices);

   deviceList = malloc(numDevices * sizeof(cl_device_id));
   if ((status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL,
                                numDevices, deviceList, NULL)) != CL_SUCCESS) {
      Warning("platform[%p]: Unable to enumerate the devices: %s\n",
              platform, CLErrString(status));
      free(deviceList);
      return;
   }

   for (ii = 0; ii < numDevices; ii++) {
      PrintDevice(deviceList[ii]);
   }

   free(deviceList);
}


int
main(int argc, char * argv[])
{
    Opts opts = { 0 };
    cl_int status;
    cl_platform_id *platformList;
    cl_uint numPlatforms;
    int ii;

    ParseOpts(&opts, argc, argv);


    if ((status = clGetPlatformIDs(0, NULL, &numPlatforms)) != CL_SUCCESS) {
       Warning("Unable to query the number of platforms: %s\n",
               CLErrString(status));
       exit(1);
    }
    printf("Found %d platform(s).\n", numPlatforms);

    platformList = malloc(sizeof(cl_platform_id) * numPlatforms);
    if ((status = clGetPlatformIDs(numPlatforms, platformList, NULL)) != CL_SUCCESS) {
       Warning("Unable to enumerate the platforms: %s\n",
               CLErrString(status));
       exit(1);
    }

    for (ii = 0; ii < numPlatforms; ii++) {
       PrintPlatform(platformList[ii]);
    }

    free(platformList);
    exit(0);
}

