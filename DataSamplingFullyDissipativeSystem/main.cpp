#include "Array2D.h"
#include "lapack.h"
#include "lapacke.h"
#include "real.h"
#include <algorithm> // std::sort
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
using std::isinf;
using std::isnan;
using std::max;
using std::vector;
using TEST::Array2D;
using TEST::real;
int nSpecies = 100;
int nReactions = 100;
double deltatime = 1;
double patankarCriticalValue = 1e-13;
int NTestedPerSample = 10;
bool conservative = true;
int NTested = 1000000;
int maximumReactants = 10;
// Coef is the coefficient for the rk scheme
int patankar(double coef, const real &density_, const real &internale,
             vector<real> &scale, const real &reaction_rate,
             const vector<vector<int>> &inputComponent,
             const vector<vector<int>> &outputComponent,
             Array2D<real> &stoichiometryCoeff, vector<real> &out) {
  double max_c_res;
  auto reactionRate = &reaction_rate;
  auto density = &density_;
  static int total_count = 0;
  static int total_nit = 0;
  static int count = 0;

  int maxReactant = 0;
  std::vector<bool> isreactant(nSpecies + 1, false);
  std::vector<real> rr(nReactions, 0);
  std::vector<int> speciesMaxPower(nSpecies + 1, 1);
  for (int re = 0; re < nReactions; ++re) {
    maxReactant = maxReactant > inputComponent[re].size()
                      ? maxReactant
                      : inputComponent[re].size();
    for (int ii = 0; ii < inputComponent[re].size(); ++ii) {
      int idx = inputComponent[re][ii];
      isreactant[idx] = true;
      speciesMaxPower[idx] =
          std::max(speciesMaxPower[idx], (int)inputComponent[re].size());
    }
  }
  real invMaxReactant = 1.0 / maxReactant;
  for (int re = 0; re < nReactions; ++re) {
    real denorm = 1.0;
    // If reaction Rate != 0, then all the scale should not be zero;
    // If reaction Rate == 0, then rr[re] = 0;
    if (reactionRate[re] != 0)
      for (int ii = 0; ii < inputComponent[re].size(); ++ii) {
        denorm *=
            pow(scale[inputComponent[re][ii]], 1.0 / inputComponent[re].size());
      }
    // If the scale is 0, the reactionRate must be 0,
    if (denorm == 0)
      rr[re] = 0;
    else
      rr[re] = reactionRate[re] * coef * deltatime / denorm;
  }
  std::vector<real> endRestultTemp(nSpecies + 1, 0);
  ++count;

  // nSpeceis + 1 equations
  // The original species
  std::vector<real> ff0(nSpecies + 1, 0);
  // Used for the equation value
  std::vector<real> ff(nSpecies + 1, 0);
  std::vector<real> reactionRateTemp(nReactions, 0);
  std::vector<real> dc(nSpecies + 1, 0);
  Array2D<double> Jacob(nSpecies + 1, nSpecies + 1);
  int *ipiv = new int[nSpecies + 1];
  double epsilon = __DBL_MIN__;

  // Firstly, generate the preliminary concentration
  // Very roughly
  double scalee = 1.0;
  for (int re = 0; re < nReactions; ++re) {
    reactionRateTemp[re] = coef * reactionRate[re];
  }
  // Initialization
  real outmax = 0;
  for (int ss = 0; ss < nSpecies; ++ss) {
    ff0[ss] = ff[ss] = density[ss];
    outmax += density[ss];
  }
  ff0[nSpecies] = ff[nSpecies] = internale;
  outmax += internale;

  std::vector<real> equationNormalizeValue(nSpecies + 1, 1.0);
  for (int eq = 0; eq < nSpecies + 1; ++eq) {
    equationNormalizeValue[eq] = std::max(ff0[eq], __DBL_MIN__);
    for (int re = 0; re < nReactions; ++re) {
      equationNormalizeValue[eq] =
          std::max(equationNormalizeValue[eq],
                   rr[re] * fabs(stoichiometryCoeff(re, eq)));
    }
  }
  std::vector<real> outMax(nSpecies + 1, outmax);

  for (int eq = 0; eq < nSpecies + 1; ++eq) {
    ff[eq] /= equationNormalizeValue[eq];
    ff0[eq] /= equationNormalizeValue[eq];
    out[eq] /= equationNormalizeValue[eq];
    scale[eq] /= equationNormalizeValue[eq];
    outMax[eq] /= equationNormalizeValue[eq];
  }
  Array2D<real> stoichiometryCoeffTemp(stoichiometryCoeff);

  // StoichiometryCoeff should also scaled.
  for (int re = 0; re < nReactions; ++re) {
    for (int eq = 0; eq < nSpecies + 1; ++eq) {
      stoichiometryCoeffTemp(re, eq) /= equationNormalizeValue[eq];
    }
  }

  for (int ss = 0; ss < nSpecies + 1; ++ss) {
    out[ss] = pow(out[ss], 1.0 / speciesMaxPower[ss]);
  }

  // Newton iteration
  std::vector<real> multiplierrr(nReactions, 0);
  std::vector<real> powerrr(nReactions, 0);
  for (int re = 0; re < nReactions; ++re) {
    powerrr[re] = 1.0 / inputComponent[re].size();
  }

  int max_nit = 10000;
  double max_c_res_0 = __DBL_MAX__;
  double theta = 1.0;
  double originaltheta = 1.0;
  std::vector<real> out0(out);
  Array2D<real> dmultiplier(nReactions, nSpecies + 1);
  dmultiplier.initialize_to_zero();
  int nit;
  real min_res = __DBL_MAX__;
  int accumulate_count = 0;
  int accumulate_count_decreasing = 0;
  int max_allowable_nondecreasing_steps = 10;
  int decreasing_step_threashold = 10;
  real factor = 0.8;
  for (nit = 0; nit < max_nit; ++nit) {
    for (int re = 0; re < nReactions; ++re) {
      multiplierrr[re] = 1.0;
      for (int ii = 0; ii < inputComponent[re].size(); ++ii) {
        int idx = inputComponent[re][ii];
        multiplierrr[re] *= pow(out[idx], speciesMaxPower[idx] * powerrr[re]);
      }
      for (int ii = 0; ii < inputComponent[re].size(); ++ii) {
        int idx = inputComponent[re][ii];
        double tmp = speciesMaxPower[idx] * powerrr[re];
        if (tmp > 1.0)
          dmultiplier(re, idx) = tmp * pow(out[idx], tmp - 1.0);
        else
          dmultiplier(re, idx) = 1.0;
        for (int jj = 0; jj < inputComponent[re].size(); ++jj) {
          if (jj == ii)
            continue;
          int idxj = inputComponent[re][jj];
          dmultiplier(re, idx) *=
              pow(out[idxj], speciesMaxPower[idxj] * powerrr[re]);
        }
      }
    }

    for (int ss = 0; ss < nSpecies + 1; ++ss) {
      ff[ss] = ff0[ss] - pow(out[ss], speciesMaxPower[ss]);
      for (int re = 0; re < nReactions; ++re) {
        ff[ss] += rr[re] * multiplierrr[re] * stoichiometryCoeffTemp(re, ss);
      }
    }

    max_c_res = -__DBL_MAX__;
    for (int i = 0; i < nSpecies + 1; ++i) {
      endRestultTemp[i] = ff[i] + pow(out[i], speciesMaxPower[i]);
      max_c_res = std::max(max_c_res, fabs(ff[i]));
    }
    // If the residual is decreasing
    // Fine
    if (max_c_res < min_res) {
      accumulate_count = 0;
      ++accumulate_count_decreasing;
      min_res = max_c_res;
      // If the residual is increasing
      // Accumulate the count
    } else {
      ++accumulate_count;
      accumulate_count_decreasing = 0;
    }
    // If the nondecreasing steps is greater than the preset threshold,
    // Decease the step size factor, i.e., theta
    if (accumulate_count > max_allowable_nondecreasing_steps) {
      accumulate_count = 0;
      min_res = max_c_res;
      theta *= factor;
    }
    // If the nonincreasing steps is greater than the preset threshold,
    // Enlarge the step size factor, i.e., theta
    if (accumulate_count_decreasing > decreasing_step_threashold) {
      theta /= factor;
      theta = std::min(theta, 1.0);
    }
    // std::cout << "Nit " << nit << " Max res " << max_c_res << std::endl;
    if (isnan(max_c_res) || isinf(max_c_res)) {
      std::cout << "NAN or INF Caught!\n";
      exit(-1);
    }
    //std::cout << "Nit is " << nit << " res is " << max_c_res << std::endl;
    if (max_c_res < patankarCriticalValue) {
      // Recalculate out[i] to keep the conservativeness;
      // It seems a littble bit difficult to keep the numerical positivity
      // rigrously although the residual is fully converged to little than
      // patankarCriticalValue; Thus still a fabs() function is applied, it
      // affects the consertivity very little.
      for (int i = 0; i < nSpecies + 1; ++i) {
        out[i] = pow(fabs(endRestultTemp[i]), 1.0 / speciesMaxPower[i]);
      }
      break;
    }
    // Then the jacobian maxtrix
    for (int ssi = 0; ssi < nSpecies + 1; ++ssi) {
      for (int ssj = 0; ssj < nSpecies + 1; ++ssj) {
        Jacob(ssi, ssj) = 0;
      }
      if (speciesMaxPower[ssi] == 1)
        Jacob(ssi, ssi) = -1.0;
      else
        Jacob(ssi, ssi) =
            -speciesMaxPower[ssi] * pow(out[ssi], speciesMaxPower[ssi] - 1);
    }
    for (int re = 0; re < nReactions; ++re) {
      for (int ii = 0; ii < inputComponent[re].size(); ++ii) {
        const int &idx = inputComponent[re][ii];
        for (int ss = 0; ss < nSpecies + 1; ++ss) {
          Jacob(ss, idx) +=
              stoichiometryCoeffTemp(re, ss) * rr[re] * dmultiplier(re, idx);
        }
      }
    }
    // Solve u^{n+1} = u^n - Jacob^-1*u^n
    lapack_int mm, info, ldu, ldvt, nn, lda;
    mm = nSpecies + 1;
    nn = nSpecies + 1;
    lda = nSpecies + 1;
    ldu = nSpecies + 1;
    ldvt = nSpecies + 1;
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, mm, nn, Jacob._data, lda, ipiv);
    info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, mm, Jacob._data, lda, ipiv);

    real min_step = 1.0;
    for (int i = 0; i < nSpecies + 1; ++i) {
      dc[i] = 0;
      for (int j = 0; j < nSpecies + 1; ++j) {
        dc[i] -= Jacob(i, j) * ff[j];
      }
      if (out[i] + dc[i] <= 0) {
        min_step = std::min(min_step, -out[i] / dc[i]);
      }
    }

    std::vector<bool> reinitialize(nSpecies + 1, false);
    bool trigger = false;
    for (int i = 0; i < nSpecies + 1; ++i) {
      // If this is not a reactant
      if (!isreactant[i]) {
        out[i] = pow(ff[i] + pow(out[i], speciesMaxPower[i]),
                     1.0 / speciesMaxPower[i]);
        // Else newton step forward;
      } else {
        out[i] += theta * dc[i];
        // Since this is a dissipative system, there is a upper bound for each
        // specie
        out[i] = std::min(out[i], pow(outMax[i], 1.0 / speciesMaxPower[i]));
        if (out[i] <= 0) {
          out[i] = 0;
          reinitialize[i] = true;
          trigger = true;
        }
      }
    }
    // If there is a negative density (although cutoff to be 0)
    // Reinitialize the cut-offed value;
    if (trigger) {
      for (int re = 0; re < nReactions; ++re) {
        multiplierrr[re] = 1.0;
        for (int ii = 0; ii < inputComponent[re].size(); ++ii) {
          int idx = inputComponent[re][ii];
          multiplierrr[re] *= pow(out[idx], speciesMaxPower[idx] * powerrr[re]);
        }
      }
      for (int ss = 0; ss < nSpecies + 1; ++ss) {
        if (!reinitialize[ss])
          continue;
        out[ss] = ff0[ss];
        for (int re = 0; re < nReactions; ++re) {
          out[ss] += rr[re] * multiplierrr[re] * stoichiometryCoeffTemp(re, ss);
        }
        out[ss] = std::min(out[ss], outMax[ss]);
        out[ss] = pow(std::max(out[ss], 0.0), 1.0 / speciesMaxPower[ss]);
      }
    }
  } // End for newton iteration
  for (int ii = 0; ii < nSpecies + 1; ++ii) {
    out[ii] = pow(out[ii], speciesMaxPower[ii]);
  }
  for (int eq = 0; eq < nSpecies + 1; ++eq) {
    scale[eq] *= equationNormalizeValue[eq];
    out[eq] *= equationNormalizeValue[eq];
  }
  delete[] ipiv;
  for (int ii = 0; ii < nSpecies + 1; ++ii) {
    if (std::isnan(out[ii]) || out[1] > 100 || out[ii] < 0) {
      std::cout << "Caught NAN or negative of patankar scheme! Not converged!"
                << std::endl;
      std::cout << "Original value is ";
      for (int ii = 0; ii < nSpecies + 1; ++ii) {
        std::cout << ff0[ii] << " ";
      }
      std::cout << "Patankar value is ";
      for (int ii = 0; ii < nSpecies + 1; ++ii) {
        std::cout << out[ii] << " ";
      }
      std::cout << std::endl;
      std::cout << "Residual for each specie is ";
      for (int ii = 0; ii < nSpecies + 1; ++ii) {
        std::cout << ff[ii] << " ";
      }
      std::cout << std::endl;
      return -1;
    }
  }
  // std::cout << " Max res " << max_c_res << std::endl;
  if (fabs(max_c_res) > patankarCriticalValue)
    return -1;
  else
    return 0;
}

double sampleTest(int argc, char **argv) {
  real density[nSpecies + 1];
  real internale;
  real reactionrate[nReactions];
  std::vector<real> scale(nSpecies + 1, 1);
  std::vector<real> out(nSpecies + 1);

  std::vector<vector<int>> inputComponent(nReactions);
  std::vector<vector<int>> outputComponent(nReactions);
  Array2D<real> stoi(nReactions, nSpecies + 1);
  std::vector<int> speciesidx(nSpecies + 1, 0);
  for (int ii = 0; ii < nSpecies + 1; ++ii) {
    speciesidx[ii] = ii;
  }
  // init the initial matrix
  for (int rr = 0; rr < nReactions; ++rr) {
    std::random_shuffle(speciesidx.begin(), speciesidx.begin() + nSpecies);
    int ninvolved = rand() % (nSpecies) + 1;
    int nreactant = rand() % ninvolved + 1;
    nreactant = std::min(nreactant, maximumReactants);
    for (int ii = 0; ii < nreactant; ++ii) {
      inputComponent[rr].push_back(speciesidx[ii]);
    }
    std::sort(inputComponent[rr].begin(), inputComponent[rr].end());
    for (int ii = nreactant; ii < ninvolved; ++ii) {
      outputComponent[rr].push_back(speciesidx[ii]);
    }
    std::sort(outputComponent[rr].begin(), outputComponent[rr].end());
    int exo = rand() % 3;
    if (exo == 0) {
      inputComponent[rr].push_back(nSpecies);
    } else if (exo == 1) {
      outputComponent[rr].push_back(nSpecies);
    }
  }
  stoi.initialize_to_zero();
  for (int rr = 0; rr < nReactions; ++rr) {
    reactionrate[rr] = 1.0;
    double total_reactant = 0;
    for (int ii = 0; ii < inputComponent[rr].size(); ++ii) {
      int idx = inputComponent[rr][ii];
      double tmp = rand() / double(RAND_MAX);
      stoi(rr, idx) = -tmp;
      total_reactant += tmp;
    }
    double total_production = 0;
    for (int ii = 0; ii < outputComponent[rr].size(); ++ii) {
      int idx = outputComponent[rr][ii];
      double tmp = rand() / double(RAND_MAX);
      stoi(rr, idx) = tmp;
      total_production += tmp;
    }
    double fraction = rand() / double(RAND_MAX);
    if(conservative) fraction = 1.0;
    for (int ii = 0; ii < outputComponent[rr].size(); ++ii) {
      int idx = outputComponent[rr][ii];
      stoi(rr, idx) *= total_reactant / total_production * fraction;
    }
  }
  for (int ss = 0; ss < nSpecies + 1; ++ss) {
    scale[ss] = rand() / double(RAND_MAX);
    density[ss] = rand() / double(RAND_MAX);
  }
  internale = density[nSpecies];
  // Read in the unconverged value;
  if (0) {
    std::ifstream fin;
    int count = 0;
    fin.open("out.dat", std::ios::binary);
    // If something error happened, the coefficients are stored in the last
    // section of the file Use while to read that section Should be improved
    // later
    if (1) {
      ++count;
      fin.read((char *)&nReactions, sizeof(nReactions));
      fin.read((char *)&nSpecies, sizeof(nSpecies));
      stoi.initialize_to_zero();
      for (int re = 0; re < nReactions; ++re) {
        int tmp;
        fin.read((char *)&tmp, sizeof(tmp));
        inputComponent[re].resize(tmp);
        for (int ii = 0; ii < inputComponent[re].size(); ++ii) {
          int tmp;
          fin.read((char *)&tmp, sizeof(tmp));
          inputComponent[re][ii] = tmp;
        }
        fin.read((char *)&tmp, sizeof(tmp));
        outputComponent[re].resize(tmp);
        for (int ii = 0; ii < outputComponent[re].size(); ++ii) {
          int tmp;
          fin.read((char *)&tmp, sizeof(tmp));
          outputComponent[re][ii] = tmp;
        }
        for (int ii = 0; ii < nSpecies + 1; ++ii) {
          double tmp;
          fin.read((char *)&tmp, sizeof(tmp));
          stoi(re, ii) = tmp;
          bool found = false;
          for (int ss = 0; ss < outputComponent[re].size(); ++ss) {
            if (outputComponent[re][ss] == ii) {
              found = true;
              break;
            }
          }
          for (int ss = 0; ss < inputComponent[re].size(); ++ss) {
            if (inputComponent[re][ss] == ii) {
              found = true;
              break;
            }
          }
          if (!found)
            stoi(re, ii) = 0;
        }
      }
      for (int ii = 0; ii < nSpecies + 1; ++ii) {
        double tmp = density[ii];
        fin.read((char *)&tmp, sizeof(tmp));
        density[ii] = tmp;
        fin.read((char *)&tmp, sizeof(tmp));
        scale[ii] = tmp;
      }
      internale = density[nSpecies];
    }
    fin.close();
  }

  int cc = 0;
  std::vector<real> out0;
  int ret;
  double totalden = 0;
  for (int ss = 0; ss < nSpecies; ++ss) {
    totalden += density[ss];
  }
  totalden += internale;
  bool error_flag = false;

  for (int ii = 0; ii < NTestedPerSample; ++ii) {
    out.resize(nSpecies + 1);
    for (int ss = 0; ss < nSpecies; ++ss) {
      out[ss] = totalden * (rand() / (RAND_MAX + 0.0));
      // out[ss] = density[ss];
    }
    out[nSpecies] = totalden * (rand() / (RAND_MAX + 0.0));
    // out[nSpecies] = internale;
    ret = patankar(1.0, density[0], internale, scale, reactionrate[0],
                   inputComponent, outputComponent, stoi, out);
    if (ret) {
      std::cout << "Fatal error found! Not converged in " << __FILE__ << ":"
                << __LINE__ << std::endl;
      ii -= 1;
      continue;
    }
    if (!ii) {
      out0 = out;
      for (int ii = 0; ii < nSpecies + 1; ++ii) {
        std::cout << std::scientific << std::setprecision(10) << out[ii] << " ";
      }
      std::cout << std::endl;
      std::cout << "Ret is " << ret << std::endl;
    } else {
      real error = 0;
      for (int ss = 0; ss < nSpecies; ++ss) {
        error += pow(out[ss] - out0[ss], 2) / (totalden * totalden);
      }
      error += pow(out[nSpecies] - out0[nSpecies], 2) /
               std::max(internale * internale, out[nSpecies] * out[nSpecies]);
      error /= (nSpecies + 1);
      error = sqrt(error);
      for (int ii = 0; ii < nSpecies + 1; ++ii) {
        std::cout << std::scientific << std::setprecision(10) << out[ii] << " ";
      }
      std::cout << std::endl;
      if (error > 1e-6) {
        std::cout << "Error happends, two solutions found" << std::endl;
        std::cout << "Ret is " << ret << std::endl;
        error_flag = true;
        goto stored;
      }
    }

  stored:
    std::ofstream fout;
    fout.open("NSpecies50NReactions50FullyDissipative.dat",
              std::ios::binary);
    fout.write((char *)&nReactions, sizeof(nReactions));
    fout.write((char *)&nSpecies, sizeof(nSpecies));
    for (int re = 0; re < nReactions; ++re) {
      int tmp = inputComponent[re].size();
      fout.write((char *)&tmp, sizeof(tmp));
      for (int ii = 0; ii < inputComponent[re].size(); ++ii) {
        int tmp = inputComponent[re][ii];
        fout.write((char *)&tmp, sizeof(tmp));
      }
      tmp = outputComponent[re].size();
      fout.write((char *)&tmp, sizeof(tmp));
      for (int ii = 0; ii < outputComponent[re].size(); ++ii) {
        int tmp = outputComponent[re][ii];
        fout.write((char *)&tmp, sizeof(tmp));
      }
      for (int ii = 0; ii < nSpecies + 1; ++ii) {
        double tmp = stoi(re, ii);
        fout.write((char *)&tmp, sizeof(tmp));
      }
    }
    for (int ii = 0; ii < nSpecies + 1; ++ii) {
      double tmp = density[ii];
      fout.write((char *)&tmp, sizeof(tmp));
      tmp = scale[ii];
      fout.write((char *)&tmp, sizeof(tmp));
      // std::cout << density[ii] << "--" << out[ii] << "\n";
    }
    fout.close();
    if (error_flag) {
      exit(-1);
    }
  }
  return ret;
}

int main(int argc, char **argv) {
  // Test sample
  // Initialize the random seed;
  srand(time(NULL));
  for (int ii = 0; ii < NTested; ++ii) {
    int ret = sampleTest(argc, argv);
    std::cout << "Testing sample " << ii << " returned value is " << ret
              << std::endl;
    if (ret != 0) {
      std::cout << "Error for the sample\n";
      exit(-1);
    }
  }
}