#!python
import install
import os

release=False
AddOption('--release', dest='release', action="store_true", default=False)
env = Environment(ENV = os.environ, RELEASE=GetOption('release'))
install.TOOL_INSTALL(env)
opts = Variables()
opts.Add(PathVariable('PREFIX', 'Directory to install under', '/usr/local/', PathVariable.PathIsDir))
opts.Add(PathVariable('JINCLUDE', 'Directory which contains jellyfish includes', None, PathVariable.PathIsDir))
opts.Add(PathVariable('JLIB', 'Directory which contains jellyfish libs', None, PathVariable.PathIsDir))
opts.Update(env)
Help(opts.GenerateHelpText(env))

# Here are our installation paths:
idir_prefix = '$PREFIX'
idir_lib    = '$PREFIX/lib'
idir_bin    = '$PREFIX/bin'
idir_etc    = '$PREFIX/etc'
idir_data   = '$PREFIX/share'
jinclude 	= '$JINCLUDE'
jlib		= '$JLIB'
Export('env idir_prefix idir_data idir_lib idir_bin idir_etc jinclude jlib')
if(env['RELEASE']==False):
	SConscript('src/SConscript', exports=['env', 'opts'])


ibin = env.InstallFiles(target=idir_bin, source='bin')
lib = env.InstallFiles(target=idir_lib, source='lib')
etc = env.InstallFiles(target=idir_etc, source='etc')
share = env.InstallFiles(target=idir_data, source='share')

#env.Alias('install', 'build') #Just installation for release with binaries
env.Alias('install', [lib, etc, share, ibin])



env.Default('build')
