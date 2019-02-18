# How to make a Kmasker Release Zip for Linux 


## Dependencies

- Computer with internet access and POSIX compatible OS (tested with MacOS Darwin 17.7.0)
- wget
- Docker Installation
- Git

### How to make it?

```bash
git clone git@github.com:tschmutzer/kmasker.git kmasker
cd kmasker
git checkout make_release_master
git submodule init
git submodule update
./make_release.sh
```
### What will happen?
Git downloads the repository and the submodules (e.g. gffread). 
The scripts downloads HolyBuildBox and starts to compile all components in src and external_src.
The compiled applications will be portable across different Linux distributions. 
Afterwards the whole kmasker folder will be zipped. 
In the end the intermediate docker containers and images are removed.

### How to update this branch before building a version?
Of course you want to build a new release with the newest commits. 
For this reason you should use 
```bash
git merge origin/master
```
before executing the make_release script. 
