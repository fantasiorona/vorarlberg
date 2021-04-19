#!python
import sys

# Create the environment and create a Compilation Database for use in VSCodium
env = DefaultEnvironment(tools=['default', 'compilation_db'])
env.CompilationDatabase()

if sys.platform.startswith('linux'):
    env.Append(CCFLAGS=['-g', '-O3'])
    program_ending = '.out'
elif sys.platform == 'darwin':
    env.Append(CCFLAGS='-g')
    program_ending = '.out'
elif sys.platform == 'win32' or sys.platform == 'msys':
    env.Append(CCFLAGS='/EHsc')
    program_ending = '.exe'
else:
    host_platform = "Unknown platform: " + sys.platform
    exit(-1)

ms_files = ["src/ms.cpp"]
tpp_files = ["src/tpp.cpp"]

Program('ms' + program_ending, ms_files)
Program('tpp' + program_ending, tpp_files)
