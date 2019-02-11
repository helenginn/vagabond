\page install Installation instructions

Vagabond uses the Meson build system, which means you will need to install meson
(requires Python 3.7 and ninja). This can be done without sudo privileges if
necessary, but it is easiest to install meson through your favourite package
manager. Ninja will come with this as a dependency.

If you need to do anything special for your system or I have a mistake
somewhere, please let me know. Vagabond is, concerningly, a little unstable on
Mac OS X - but the command line version usually works.

## Installing dependencies

The exact packages may vary between yum/dnf/apt-get and I don't know the
commands on all the systems yet! These commands will need sudo privileges.

One dependency, `libpng` is probably already available on your system and
does not need to be installed with your package manager. 

Once you have meson, please install the Fast Fourier Transform library
development version ("-devel" versions allow pkg-config to find the include
headers and libraries for compiling other programs from scratch):

```
dnf install fftw-devel
```

Not sure for apt-get, but for Brew (Mac OS X) (needs checking):

```
brew install fftw
```

If you want to be able to run the GUI, you will need Qt5.

```
dnf install qt5-devel
```

If you are using apt-get, I think you will need this instead:

```
apt-get install qt3d5-dev
```

Using brew,

```
brew install qt5
```

## Installing Vagabond

Firstly, clone the repository and enter the directory

```
git clone https://www.github.com/helenginn/vagabond
cd vagabond
```

On Linux, we can configure a build directory and compile it:

```
meson --buildtype=release build/current
```

At this point, you will be notified of any missing libraries which meson
could not find. It will also only compile the command line version if it
was unable to find Qt5. 

Then, compile the libraries and binaries:

```
ninja -C build/current
```

On Mac OS X, we need to pass a few flags to the compiler, unfortunately.

```
CXXFLAGS='-mmacosx-version-min=10.7 --std=c++0x -stdlib=libc++' meson build/current
```
