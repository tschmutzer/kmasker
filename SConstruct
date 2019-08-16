#!python
import install

env = Environment(ENV = os.environ)
install.TOOL_INSTALL(env)
opts = Variables()
opts.Add(PathVariable('PREFIX', 'Directory to install under', '/usr', PathVariable.PathIsDir))
opts.Add(PathVariable('JINCLUDE', 'Directory which contains jellyfish includes', None, PathVariable.PathIsDir))
opts.Update(env)
Help(opts.GenerateHelpText(env))

# Here are our installation paths:
idir_prefix = '$PREFIX'
idir_lib    = '$PREFIX/lib'
idir_bin    = '$PREFIX/bin'
idir_etc    = '$PREFIX/etc'
idir_data   = '$PREFIX/share'
jinclude 	= '$JINCLUDE'
Export('env idir_prefix idir_data idir_lib idir_bin idir_etc jinclude')

SConscript('src/SConscript', exports=['env', 'opts'])
bin = env.InstallFiles(target=idir_bin, source='bin')
lib = env.InstallFiles(target=idir_lib, source='lib')
etc = env.InstallFiles(target=idir_etc, source='etc')
share = env.InstallFiles(target=idir_data, source='share')

env.Alias('install', 'build')
env.Alias('install', lib)
env.Alias('install', etc)
env.Alias('install', share)
env.Alias('install', bin)

env.Default('build')
