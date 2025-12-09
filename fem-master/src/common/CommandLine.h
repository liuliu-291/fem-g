// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_COMMANDLINE
#define H_GMSHFEM_COMMANDLINE

#include <functional>
#include <string>
#include <vector>

namespace gmshfem::common
{
  class ArgsManager;
}

namespace gmshfem::common
{


  namespace CommandLine
  {

    enum class Type {
      Search,
      Help,
      Unary,
      Binary,
      Ternary
    };

    struct Parameter {
      enum class Type {
        Real,
        Integer,
        String
      };

      double d;
      int i;
      std::string s;

      Parameter() :
        d(0.), i(0), s() {}
      ~Parameter() {}
      Parameter &operator=(const Parameter &other)
      {
        d = other.d;
        i = other.i;
        s = other.s;

        return *this;
      }
    };


  } // namespace CommandLine

  class CommandLineInterface
  {
   protected:
    const std::string _name;
    const std::string _help;
    const bool _kill;
    const bool _pre;
    const std::vector< CommandLine::Parameter::Type > _parametersType;
    std::vector< CommandLine::Parameter > _parameters;

   public:
    CommandLineInterface(const std::string &name, const std::string &help, const bool kill, const bool pre, const std::vector< CommandLine::Parameter::Type > &parametersType);
    virtual ~CommandLineInterface();

    const std::string &name() const;
    const std::string &help() const;

    bool kill() const;
    bool pre() const;
    bool post() const;

    virtual CommandLine::Type type() const = 0;

    virtual CommandLine::Parameter::Type
    getParameterType(const unsigned int i = 0);
    virtual std::string getParameterTypeName(const unsigned int i = 0);
    virtual void setParameter(const CommandLine::Parameter &parameter, const unsigned int i = 0);

    virtual void run(ArgsManager *argsManager) const = 0;

    unsigned int nbrParameters() const;

    virtual CommandLineInterface *copy() const = 0;
  };

  class SearchCL final : public CommandLineInterface
  {
   public:
    explicit SearchCL(const std::string &name);
    ~SearchCL();

    virtual CommandLine::Type type() const override;

    virtual void run(ArgsManager *argsManager) const override;

    CommandLineInterface *copy() const override;
  };

  class UnaryCL : public CommandLineInterface
  {
   private:
    void (*_action)(ArgsManager *argsManager);

   public:
    UnaryCL(const std::string &name, const std::string &help, const bool kill, const bool pre, void (*action)(ArgsManager *argsManager));
    ~UnaryCL();

    virtual CommandLine::Type type() const override;

    virtual void run(ArgsManager *argsManager) const override;

    CommandLineInterface *copy() const override;

    void action(void (*action)(ArgsManager *argsManager));
  };

  class UnaryHelpCL final : public UnaryCL
  {
   public:
    UnaryHelpCL();
    ~UnaryHelpCL();

    void run(ArgsManager *argsManager) const override;

    CommandLine::Type type() const override;

    CommandLineInterface *copy() const override;
  };

  class BinaryCL final : public CommandLineInterface
  {
   private:
    void (*_action)(ArgsManager *argsManager, const CommandLine::Parameter &parameter);

   public:
    BinaryCL(const std::string &name, const std::string &help, const bool kill, const bool pre, const CommandLine::Parameter::Type &ParameterType, void (*action)(ArgsManager *argsManager, const CommandLine::Parameter &parameter));
    ~BinaryCL();

    virtual CommandLine::Type type() const override;

    virtual void run(ArgsManager *argsManager) const override;

    CommandLineInterface *copy() const override;

    void action(void (*action)(ArgsManager *argsManager, const CommandLine::Parameter &parameter));
  };


} // namespace gmshfem::common

// std::less

namespace std
{


  template<>
  struct less< gmshfem::common::CommandLineInterface * > {
    bool operator()(const gmshfem::common::CommandLineInterface *a, const gmshfem::common::CommandLineInterface *b) const
    {
      return a->name() < b->name();
    }
  };


} // namespace std

#endif // H_GMSHFEM_COMMANDLINE
