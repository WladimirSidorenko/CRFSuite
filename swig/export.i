%module(directors="1") crfsuite

%include <std_string.i>
%include <std_vector.i>
%include <exception.i>

%{
#include "crfsuite_api.hpp"
%}

%include "crfsuite_api.hpp"

namespace CRFSuite {
  %template(Item) ::std::vector<Attribute>;
   %template(ItemSequence) ::std::vector<Item>;
   %template(StringList) ::std::vector<std::string>;
 }
typedef CRFSuite::StringList StringList;

%feature("director") Trainer;

%exception {
    try {
        $action
    } catch(const std::invalid_argument& e) {
        SWIG_exception(SWIG_IOError, e.what());
    } catch(const std::runtime_error& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch(...) {
        SWIG_exception(SWIG_RuntimeError,"Unknown exception");
    }
}
