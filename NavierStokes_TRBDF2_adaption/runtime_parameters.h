// We start by including all the necessary deal.II header files
//
#include <deal.II/base/parameter_handler.h>

// @sect{Run time parameters}
//
// Since our method has several parameters that can be fine-tuned we put them
// into an external file, so that they can be determined at run-time.
//
namespace RunTimeParameters {
  using namespace dealii;

  class Data_Storage {
  public:
    Data_Storage();

    void read_data(const std::string& filename);

    double initial_time;
    double final_time;

    double Reynolds;
    double Cs2;
    double dt;

    unsigned int n_refines;            /*--- Number of refinements ---*/
    unsigned int max_loc_refinements; /*--- Number of maximum local refinements allowed ---*/
    unsigned int min_loc_refinements; /*--- Number of minimum local refinements allowed
                                            once reached that level ---*/
    bool big_domain;                  /*--- Flag which domain one wants consider ---*/
    bool square_cylinder;             /*--- Flag to determine whether one wants
                                            to simulate the square cylinder or not ---*/
    bool empty;                       /*--- Flag to define an empty channel ---*/

    /*--- Parameters related to the linear solver ---*/
    unsigned int max_iterations;
    double       eps;
    double       tolerance_fixed_point;

    bool         verbose;
    unsigned int output_interval;

    std::string dir; /*--- Auxiliary string variable for output storage ---*/

    unsigned int refinement_iterations; /*--- Auxiliary variable about how many steps perform remeshing ---*/

    /*--- Auxiliary parameters related to restart ---*/
    bool         restart;
    bool         save_for_restart;
    unsigned int step_restart;
    double       time_restart;
    bool         as_initial_conditions;
    bool         modify_Reynolds;

  protected:
    ParameterHandler prm;
  };

  // In the constructor of this class we declare all the parameters in suitable (but arbitrary) subsections.
  //
  Data_Storage::Data_Storage(): initial_time(0.0),
                                final_time(1.0),
                                Reynolds(1.0),
                                dt(5e-4),
                                n_refines(0),
                                max_loc_refinements(0),
                                min_loc_refinements(0),
                                big_domain(false),
                                square_cylinder(false),
                                empty(false),
                                max_iterations(1000),
                                eps(1e-12),
                                tolerance_fixed_point(1e-6),
                                verbose(true),
                                output_interval(15),
                                refinement_iterations(0),
                                restart(false),
                                save_for_restart(false),
                                step_restart(0),
                                time_restart(0.0),
                                as_initial_conditions(false),
                                modify_Reynolds(false) {
    prm.enter_subsection("Physical data");
    {
      prm.declare_entry("initial_time",
                        "0.0",
                        Patterns::Double(0.0),
                        " The initial time of the simulation. ");
      prm.declare_entry("final_time",
                        "1.0",
                        Patterns::Double(0.0),
                        " The final time of the simulation. ");
      prm.declare_entry("Reynolds",
                        "1.0",
                        Patterns::Double(0.0),
                        " The Reynolds number. ");
      prm.declare_entry("Cs2",
                        "0.01",
                        Patterns::Double(0.0),
                        " The Smagorinsky coefficient. ");
    }
    prm.leave_subsection();

    prm.enter_subsection("Time step data");
    {
      prm.declare_entry("dt",
                        "5e-4",
                        Patterns::Double(0.0),
                        " The time step size. ");
      prm.declare_entry("time_restart",
                        "5e-4",
                        Patterns::Double(0.0),
                        " The time of restart. ");
    }
    prm.leave_subsection();

    prm.enter_subsection("Space discretization");
    {
      prm.declare_entry("n_of_refines",
                        "100",
                        Patterns::Integer(0, 1500),
                        " The number of cells we want on each direction of the mesh. ");
      prm.declare_entry("max_loc_refinements",
                        "4",
                         Patterns::Integer(0, 10),
                         " The number of maximum local refinements. ");
      prm.declare_entry("min_loc_refinements",
                        "2",
                         Patterns::Integer(0, 10),
                         " The number of minimum local refinements. ");
      prm.declare_entry("big_domain",
                        "false",
                        Patterns::Bool(),
                        " This flag decides if the domain is chosen small or big. ");
      prm.declare_entry("square_cylinder",
                        "false",
                        Patterns::Bool(),
                        " This flag decides if we consider square or round cylinders. ");
      prm.declare_entry("empty",
                        "false",
                        Patterns::Bool(),
                        " This flag decides if we consider an empty channel. ");
    }
    prm.leave_subsection();

    prm.enter_subsection("Data solve");
    {
      prm.declare_entry("max_iterations",
                        "1000",
                        Patterns::Integer(1, 30000),
                        " The maximal number of iterations linear solvers must make. ");
      prm.declare_entry("eps",
                        "1e-12",
                        Patterns::Double(0.0),
                        " The stopping criterion. ");
      prm.declare_entry("step_restart",
                        "1",
                         Patterns::Integer(1, 100000000),
                         " The step at which restart occurs");
      prm.declare_entry("tolerance_fixed_point",
                        "1e-8",
                         Patterns::Double(0.0),
                         " The tolerance of the fixed point loop");
    }
    prm.leave_subsection();

    prm.declare_entry("refinement_iterations",
                      "0",
                      Patterns::Integer(0),
                      " This number indicates how often we need to "
                      "refine the mesh");

    prm.declare_entry("saving directory", "SimTest");

    prm.declare_entry("verbose",
                      "true",
                      Patterns::Bool(),
                      " This indicates whether the output of the solution "
                      "process should be verbose. ");

    prm.declare_entry("output_interval",
                      "1",
                      Patterns::Integer(1),
                      " This indicates between how many time steps we print "
                      "the solution. ");

    prm.declare_entry("restart",
                      "false",
                      Patterns::Bool(),
                      " This indicates whether we are in presence of a "
                      "restart or not. ");
    prm.declare_entry("save_for_restart",
                      "false",
                      Patterns::Bool(),
                      " This indicates whether we want to save for possible "
                      "restart or not. ");
    prm.declare_entry("as_initial_conditions",
                      "false",
                      Patterns::Bool(),
                      " This indicates whether restart is used as initial condition "
                      "or to continue the simulation. ");
    prm.declare_entry("modify_Reynolds",
                      "false",
                      Patterns::Bool(),
                      " This indicates whether we want to change manually the "
                      " Reynolds number. ");
  }

  // We need now a routine to read all declared parameters in the constructor
  //
  void Data_Storage::read_data(const std::string& filename) {
    std::ifstream file(filename);
    AssertThrow(file, ExcFileNotOpen(filename));

    prm.parse_input(file);

    prm.enter_subsection("Physical data");
    {
      initial_time = prm.get_double("initial_time");
      final_time   = prm.get_double("final_time");
      Reynolds     = prm.get_double("Reynolds");
      Cs2          = prm.get_double("Cs2");
    }
    prm.leave_subsection();

    prm.enter_subsection("Time step data");
    {
      dt           = prm.get_double("dt");
      time_restart = prm.get_double("time_restart");
    }
    prm.leave_subsection();

    prm.enter_subsection("Space discretization");
    {
      n_refines           = prm.get_integer("n_of_refines");
      max_loc_refinements = prm.get_integer("max_loc_refinements");
      min_loc_refinements = prm.get_integer("min_loc_refinements");
      big_domain          = prm.get_bool("big_domain");
      square_cylinder     = prm.get_bool("square_cylinder");
      empty               = prm.get_bool("empty");
    }
    prm.leave_subsection();

    prm.enter_subsection("Data solve");
    {
      max_iterations        = prm.get_integer("max_iterations");
      eps                   = prm.get_double("eps");
      step_restart          = prm.get_integer("step_restart");
      tolerance_fixed_point = prm.get_double("tolerance_fixed_point");
    }
    prm.leave_subsection();

    dir = prm.get("saving directory");

    refinement_iterations = prm.get_integer("refinement_iterations");

    verbose = prm.get_bool("verbose");

    output_interval = prm.get_integer("output_interval");

    /*--- Read parameters related to restart ---*/
    restart               = prm.get_bool("restart");
    save_for_restart      = prm.get_bool("save_for_restart");
    as_initial_conditions = prm.get_bool("as_initial_conditions");
    modify_Reynolds       = prm.get_bool("modify_Reynolds");
  }

} // namespace RunTimeParameters
