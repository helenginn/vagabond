#ifndef __fuck_cov__MtzFFTPtr__
#define __fuck_cov__MtzFFTPtr__

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

class MtzFFT;
typedef boost::shared_ptr<MtzFFT> MtzFFTPtr;
#define ToMtzFFTPtr(a) (boost::static_pointer_cast<MtzFFT>((a)))


#endif
