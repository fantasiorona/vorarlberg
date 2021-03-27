#!python

# Create the environment and create a Compilation Database for use in VSCodium
env = DefaultEnvironment(tools=['default', 'compilation_db'])
env.CompilationDatabase()

env.Append(CCFLAGS='-g')

Program('ms.out', Glob('src/*.cpp'))
