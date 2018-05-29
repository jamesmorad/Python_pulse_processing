#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "PulseFinder.h"
#include "SplitGaussian.h"
using namespace boost::python;

template<class T> boost::python::list std_vector_to_py_list(const std::vector<T>& v)
{
    boost::python::object get_iter = boost::python::iterator<std::vector<T> >();
    boost::python::object iter = get_iter(v);
    boost::python::list l(iter);
    return l;
}

template< class T >
inline
std::vector< T > to_std_vector( const boost::python::object& iterable )
{
    return std::vector< T >( boost::python::stl_input_iterator< T >( iterable ),
                             boost::python::stl_input_iterator< T >( ) );
}


BOOST_PYTHON_MODULE(PulseFinder_BlackBox){
    
    //! python access to stl integer vectors
    class_< std::vector<float> >("vectorFloat")
        .def(vector_indexing_suite<std::vector<float> >())
    ;

    class_< std::vector<int> >("vectorInt")
        .def(vector_indexing_suite<std::vector<int> >())
    ;


    def("std_vector_to_py_list_float",std_vector_to_py_list<float>);
    def("std_vector_to_py_list_int",std_vector_to_py_list<int>);
    def("to_std_vector",to_std_vector<float>);



    class_<PulseFinder>("PulseFinder")
        .def("Initilize",&PulseFinder::Initilize)
        .def("Execute",&PulseFinder::Execute)
        .def_readwrite("kz_samples", &PulseFinder::kz_samples)
        .def_readwrite("kz_iterations", &PulseFinder::kz_iterations)
        .def_readwrite("min_max_ratio", &PulseFinder::min_max_ratio)
        .def_readwrite("max_threshold", &PulseFinder::max_threshold)
        .def_readwrite("max_seperation", &PulseFinder::max_seperation)
        .def_readwrite("nsigma", &PulseFinder::nsigma)
        .def_readwrite("pre_buffer_samples", &PulseFinder::pre_buffer_samples)
        .def_readwrite("post_buffer_samples", &PulseFinder::post_buffer_samples)
        .def_readwrite("se_amp", &PulseFinder::se_amp)
        .def_readwrite("se_area_mean_phe", &PulseFinder::se_area_mean_phe)
        .def_readwrite("se_area_sigma_phe", &PulseFinder::se_area_sigma_phe)
        .def_readwrite("se_width_mean_samples", &PulseFinder::se_width_mean_samples)
        .def_readwrite("se_width_sigma_samples", &PulseFinder::se_width_sigma_samples)
        .def_readwrite("s1AMean", &PulseFinder::s1AMean)
        .def_readwrite("s1ASigma", &PulseFinder::s1ASigma)
        .def_readwrite("s1AOffset", &PulseFinder::s1AOffset)
        .def_readwrite("s1HMean", &PulseFinder::s1HMean)
        .def_readwrite("s1HSigma", &PulseFinder::s1HSigma)
        .def_readwrite("s1HOffset", &PulseFinder::s1HOffset)
        .def_readwrite("s2AMean", &PulseFinder::s2AMean)
        .def_readwrite("s2ASigma", &PulseFinder::s2ASigma)
        .def_readwrite("s2AOffset", &PulseFinder::s2AOffset)
        .def_readwrite("s2HMean", &PulseFinder::s2HMean)
        .def_readwrite("s2HSigma", &PulseFinder::s2HSigma)
        .def_readwrite("s2HOffset", &PulseFinder::s2HOffset)
        .def_readwrite("s1SDscale", &PulseFinder::s1SDscale)
        .def_readwrite("s2SDscale", &PulseFinder::s2SDscale)
        .def_readwrite("wave", &PulseFinder::wave)
        .def_readwrite("smooth", &PulseFinder::smooth)
        .def_readwrite("borders", &PulseFinder::borders)
        .def_readonly("pars", &PulseFinder::pars)

    ;
  
    
}