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

target_name = 'ms' + program_ending

Program(target_name, Glob('src/*.cpp'))
