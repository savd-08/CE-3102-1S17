#ifndef _EXCEPTION
#define _EXCEPTION

#include <exception>
#include <string>

/**
* Excepcion general
*
*/
class Exception: public std::exception
{
public:
    /**
     *  Constructor con puntero de string
     */
    explicit Exception(const char* message):msg_(message){}

    /**
     *  Constructor con objeto string
     */
    explicit Exception(const std::string& message):msg_(message){}

    /**
     * Destructor.
     */
    virtual ~Exception() throw (){}

    /**
     *  Retorna puntero de la descripcion de la excepcion
     */
    virtual const char* what() const throw (){
       return msg_.c_str();
    }

protected:
    /**
     * Descricion de excepcion
     */
    std::string msg_;
};

#endif
