#include "gtest/gtest.h"
#include "AppFactory.h"
#include "MooseInit.h"
#include "Moose.h"
#include "MooseApp.h"
#include "CarbonDioxideApp.h"

PerfLog Moose::perf_log("gtest");

GTEST_API_ int
main(int argc, char ** argv)
{
  // gtest removes (only) its args from argc and argv - so this  must be before moose init
  testing::InitGoogleTest(&argc, argv);

  MooseInit init(argc, argv);
  registerApp(CarbonDioxideApp);
  Moose::_throw_on_error = true;

  return RUN_ALL_TESTS();
}
