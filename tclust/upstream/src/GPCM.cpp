/**
 * GPCM - Parsimonious constraints implementation
 *
 */

//  Rcpp::compileAttributes("C:/users/valen/onedrive/myrepo/r/tclust")

// [[Rcpp::depends(RcppArmadillo)]]

#include "tclust_types.h"

#include <Rcpp.h>

using namespace Rcpp;

arma::mat restr2Eigenv(arma::mat autovalues, arma::vec ni_ini, double factor_e, double zero_tol);

arma::cube reconstructSigma(const arma::vec& lmd,
                            const arma::cube& OMG,
                            const arma::mat& GAM)
{
    const unsigned int p = OMG.n_rows;
    const unsigned int k = OMG.n_slices;

    arma::cube Sigma(p, p, k);

    for (unsigned int j = 0; j < k; ++j)
    {
        Sigma.slice(j) =
            lmd(j) *
            OMG.slice(j) *
            arma::diagmat(GAM.col(j)) *
            OMG.slice(j).t();
    }

    return Sigma;
}

//  Finds initial common rotation matrix
//
// This procedures is called when **E that is when a common rotation matrix
// is imposed. The main purpose of this function is to find the an initial
// estimate of the common rotation matrix.
//
Rcpp::List initR(const arma::cube& SigmaB,
                 const arma::vec& niini,
                 int p,
                 int k,
                 const std::string& pars)
{
    arma::vec lmd(k);

    if (pars[0] == 'V')
    {

// DEBUG
//Rcout << "In initR()" << std::endl;
//SigmaB.print("SigmaB");

        // VT::22.07.2026 - The problem with zero and nan lmd in case of VVE - replace det() by eig_sym()
        //  arma::vec eigval;
        //  arma::mat eigvec;   

        for (int j = 0; j < k; ++j)
        {
            // VT::22.07.2026 - The problem with zero and nan lmd in case of VVE - replace det() by eig_sym()
            //  arma::eig_sym(eigval, eigvec, SigmaB.slice(j));
            //  lmd(j) = std::pow(arma::prod(eigval), 1.0/p);
            
            // DEBUG
            //double tolerance = 1.0 * arma::datum::eps; 
            //bool xcond = arma::rcond(SigmaB.slice(j)) >= tolerance;
            //Rcout << "initR: " << j << "   " << arma::rcond(SigmaB.slice(j)) << "   " << tolerance << "   " << xcond << std::endl;
            //eigval.print("eigval");
            
            // Rcout << "initR: " << j << "   " << arma::det(SigmaB.slice(j)) << "   " << std::pow(arma::det(SigmaB.slice(j)), 1.0 / p) << std::endl;
            lmd(j) = std::pow(arma::det(SigmaB.slice(j)), 1.0 / p);
            
            //if(lmd(j) == 0)
            //    lmd(j) = 1;
        }
    }
    else
    {
        lmd.ones();
    }

// DEBUG
//Rcout << "In initR()" << std::endl;
//lmd.print("lmd");

    arma::mat Sw(p, p, arma::fill::zeros);

    double n = arma::sum(niini);

    for (int j = 0; j < k; ++j)
    {
        Sw += (niini(j) / n) *
              (1.0 / lmd(j)) *
              SigmaB.slice(j);
    }

    arma::mat Omega2D;

    if (!Sw.is_finite())
    {
        Omega2D.eye(p, p);
    }
    else
    {
        arma::vec eigval;
        arma::mat eigvec;

        arma::eig_sym(eigval, eigvec, Sw);

        Omega2D = arma::fliplr(eigvec);
    }

    arma::cube Omega(p, p, k);

    for (int j = 0; j < k; ++j)
        Omega.slice(j) = Omega2D;

    return Rcpp::List::create(
        Rcpp::Named("lmd") = lmd,
        Rcpp::Named("Omega") = Omega,
        Rcpp::Named("Omega2D") = Omega2D
    );
}

//  Computes updated common rotation matrix when shapes are equal
//
//  This routine is called when the parameterization is VEE, that is when
//  equal shape and equal rotation are imposed and we have varying
//  determinants.
//
Rcpp::List cpcE(const arma::vec& lmd,
                const arma::cube& SigmaB,
                const arma::vec& niini)
{
    const int p = SigmaB.n_rows;
    const int k = SigmaB.n_slices;

    arma::mat Sigma(p, p, arma::fill::zeros);

    double n = arma::sum(niini);

    for (int j = 0; j < k; ++j)
    {
        Sigma +=
            (niini(j) / n) *
            (1.0 / lmd(j)) *
            SigmaB.slice(j);
    }

    arma::vec eigval;
    arma::mat eigvec;

    arma::eig_sym(eigval, eigvec, Sigma);

    arma::mat Omega2D = arma::fliplr(eigvec);

    arma::cube Omega(p, p, k);

    for (int j = 0; j < k; ++j)
        Omega.slice(j) = Omega2D;

    return Rcpp::List::create(
        Rcpp::Named("Omega") = Omega,
        Rcpp::Named("Omega2D") = Omega2D
    );
}

//  Computes updated common rotation matrix when shapes are different
//
//  This routine is called when the parameterization is *VE, that is when
//  variable shape is assumed but equal rotation is imposed. This routine is
//  based on the algorithm described in McNicholas and Browne (2014)
//
Rcpp::List cpcV(const arma::vec& lmd,
                const arma::mat& GAM,
                arma::mat Omega2D,
                const arma::cube& Wk,
                const arma::vec& wk,
                double tolR,
                int maxiterR)
{
    const int p = Wk.n_rows;
    const int k = Wk.n_slices;

    arma::mat OmegaOld(p, p, arma::fill::zeros);

    double diff = std::numeric_limits<double>::infinity();

    int iter = 0;

    while ((diff > tolR) && (iter < maxiterR))
    {
        ++iter;

        arma::mat F(p, p, arma::fill::zeros);

        for (int j = 0; j < k; ++j)
        {
            arma::mat Ginv = arma::diagmat(1.0 / GAM.col(j));

            double fac = std::pow(lmd(j), -1.0 / p);

            F += fac *
                 Ginv *
                 Omega2D.t() *
                 Wk.slice(j);

            F -= fac *
                 wk(j) *
                 Ginv *
                 Omega2D.t();
        }

        arma::mat U, V;
        arma::vec s;

        arma::svd(U, s, V, F);

        Omega2D = V * U.t();

        if (iter > 1)
        {
            arma::mat A = Omega2D.t() * OmegaOld;
            diff = (p - arma::accu(arma::square(A))) / p;
        }

        OmegaOld = Omega2D;
    }

    arma::cube Omega(p, p, k);

    for (int j = 0; j < k; ++j)
        Omega.slice(j) = Omega2D;

    return Rcpp::List::create(
        Rcpp::Named("Omega") = Omega,
        Rcpp::Named("Omega2D") = Omega2D
    );
}

// Computes constrained Gamma (shape) matrix
//
// The purpose is to find the new constrained shape matrix.
// This routine is called when modeltype(2) = V (that is in presence of
// varying shape matrices). In order to impose both shw inside each group
// and shb between groups, a number of iterations specified by input
// parameter maxiterS and a stopping condition given by itertol are necessary .
//
arma::mat restrshapecore(const arma::mat& GAMini,
                         const arma::vec& niini,
                         double shw,
                         double shb,
                         double zerotol,
                         int maxiterS,
                         double tolS,
                         bool sortsh = true)
{
    const unsigned int p = GAMini.n_rows;
    const unsigned int K = GAMini.n_cols;

    arma::mat lamGAMc = GAMini;

    // DEBUG
    //Rcout << "Entering restrshapecore()" << std::endl;
    //Rcout << "shw" <<shw << "shb" <<shb << std::endl;
    //GAMini.print("GAMini");
    //niini.print("niini");

    //------------------------------------------------------------------
    // Step 1 : within-group restriction (columns)
    //------------------------------------------------------------------

    for (unsigned int j = 0; j < K; ++j)
    {
        arma::vec g = GAMini.col(j);

        // if (g.max() / g.min() > shw)
        double mn = std::max(g.min(), zerotol);
        if (g.max()/mn > shw)
        {
            arma::mat tmp = restr2Eigenv(g, arma::vec{1.0}, shw, zerotol);
            lamGAMc.col(j) = tmp.col(0);
        }
    }

    // DEBUG
    //Rcout << "Step 1 : within-group restriction (columns)" << std::endl;
    //lamGAMc.print("lamGAMc");
    
    //------------------------------------------------------------------
    // Main iteration
    //------------------------------------------------------------------

    arma::mat GAM(p, K);
    arma::mat GAMsor(p, K);
    arma::mat GAMctr(p, K);
    arma::mat GAMctrSRT(p, K);

    arma::umat Ord;

    double diffGAM = std::numeric_limits<double>::infinity();

    int iter = 0;

    while ((diffGAM > tolS) && (iter < maxiterS))
    {
        ++iter;

        GAM = lamGAMc;

        arma::vec GAMold = arma::vectorise(GAM);

        //--------------------------------------------------------------
        // det(Gamma)=1
        //--------------------------------------------------------------

        arma::rowvec scale = arma::pow(arma::prod(GAM,0), 1.0/p);

        scale.replace(0.0,1.0);

        GAM.each_row() /= scale;

        GAM.replace(0.0,1.0);

        //--------------------------------------------------------------
        // sort each column
        //--------------------------------------------------------------

        if(sortsh)
        {
            // DEBUG
            //Rcout << "sortsh is TRUE - sort each column of GAM" << std::endl;
            //GAM.print("GAM");
            
            Ord.set_size(p,K);

            for(unsigned int j=0;j<K;++j)
            {
                //Ord.col(j)=arma::sort_index(GAM.col(j));
                // GAMsor.col(j)=GAM.col(j).elem(Ord.col(j));
                
                arma::vec g = GAM.col(j);
                Ord.col(j) = arma::sort_index(g);
                GAMsor.col(j) = g.elem(Ord.col(j));
            }
        }
        else
        {
            GAMsor=GAM;
        }

        //--------------------------------------------------------------
        // between-group restriction (rows)
        //--------------------------------------------------------------

        GAMctr=GAMsor;

        for(unsigned int i=0;i<p;++i)
        {
            arma::rowvec r=GAMsor.row(i);

            if(r.max()/r.min()>shb)
            {
                arma::mat tmp = restr2Eigenv(r.t(), niini, shb, zerotol);

                GAMctr.row(i)=tmp.col(0).t();
            }
        }

        //--------------------------------------------------------------
        // restore original ordering
        //--------------------------------------------------------------

        if(sortsh)
        {
            for(unsigned int j=0;j<K;++j)
            {

                // Replace te nested loop below.
                // This is the Armadillo idiom for applying the inverse permutation 
                //  and is both clearer and typically faster.

                // for(unsigned int i=0;i<p;++i)
                //    GAMctrSRT(Ord(i,j),j)=GAMctr(i,j);
                
                arma::vec tmp(p);
                tmp.elem(Ord.col(j)) = GAMctr.col(j);
                GAMctrSRT.col(j) = tmp;
                
                
            }

            lamGAMc=GAMctrSRT;
        }
        else
        {
            lamGAMc=GAMctr;
        }

        arma::vec GAMnew=arma::vectorise(lamGAMc);

        diffGAM=
            arma::accu(arma::square(GAMnew-GAMold))/
            arma::accu(arma::square(GAMold));
    }

    return lamGAMc;
}

// Produces the restricted shape matrix for the 14 GPCM
//
//
// The purpose of this routine is to produce the constrained shape matrix
// $\Gamma$.
// This routine copes with the second of the 3 letters of the model type. It
// deals with the cases in which the second letter is E, or I or V. If the
// second letter is V procedure restrshapecore is invoked and both (within
// groups) cshw, and (between groups) cshb constraints are imposed. If the
// second letter of model type is E just cshw is used. If the second letter
// is I, GAMc becomes a matrix of ones.
//
arma::mat restrshapeGPCM(const arma::vec& lmd,
                         const arma::cube& Omega,
                         const arma::cube& SigmaB,
                         const arma::vec& niini,
                         const GPCMPars& pa)
{
    const unsigned int p = pa.p;
    const unsigned int k = pa.k;

    //------------------------------------------------------------------
    // Common shape
    //------------------------------------------------------------------

    if ((pa.pars[1] == 'E') || (pa.shb == 1.0))
    {
        arma::mat GAMpooled(p,p,arma::fill::zeros);

        double n = arma::sum(niini);

        for(unsigned int j=0;j<k;++j)
        {
            arma::mat G =
                (niini(j)/n)*
                (1.0/lmd(j))*
                Omega.slice(j).t()*
                SigmaB.slice(j)*
                Omega.slice(j);

            if(!G.is_finite())
                G.zeros();

            GAMpooled += G;
        }

        arma::mat ev = GAMpooled.diag();

        arma::mat GAMpooledc = restr2Eigenv(ev, arma::vec{1.0}, pa.shw, pa.zerotol);

        double scale = std::pow(arma::prod(GAMpooledc.col(0)), 1.0/p);

        if(scale==0.0)
            scale=1.0;

        arma::vec shape = GAMpooledc.col(0)/scale;

        arma::mat GAMc(p,k);

        for(unsigned int j=0;j<k;++j)
            GAMc.col(j)=shape;

        return GAMc;
    }

    //------------------------------------------------------------------
    // Identity shape
    //------------------------------------------------------------------

    if ((pa.pars[1] == 'I') || (pa.shw == 1.0))
    {
        return arma::mat(p,k,arma::fill::ones);
    }

    //------------------------------------------------------------------
    // Variable shape
    //------------------------------------------------------------------

    arma::mat GAM(p,k);

    for(unsigned int j=0;j<k;++j)
    {
        arma::mat tmp =
            Omega.slice(j).t()*
            SigmaB.slice(j)*
            Omega.slice(j);

        GAM.col(j)=tmp.diag()/lmd(j);
    }

    return restrshapecore(GAM,
                          niini,
                          pa.shw,
                          pa.shb,
                          pa.zerotol,
                          pa.maxiterS,
                          pa.tolS,
                          pa.sortsh!=0);
}

// Applies determinant restrictions for the 14 GPCM
//
//  This routine applies the constraints to the determinants using the
//  specification contained in field pa.cdet of input structure pa.
//
//
arma::vec restrdeterGPCM(const arma::mat& GAM,
                         const arma::cube& OMG,
                         const arma::cube& SigmaB,
                         const arma::vec& niini,
                         const GPCMPars& pa)
{
    const unsigned int k = pa.k;
    const unsigned int p = pa.p;

    arma::vec lmd(k);

    //------------------------------------------------------------------
    // Unconstrained determinants
    //------------------------------------------------------------------

    for (unsigned int j = 0; j < k; ++j)
    {
        arma::mat tmp =
            OMG.slice(j).t() *
            SigmaB.slice(j) *
            OMG.slice(j);

        lmd(j) = arma::accu(tmp.diag() / GAM.col(j)) / static_cast<double>(p);
    }

    //------------------------------------------------------------------
    // Apply restriction if needed
    //------------------------------------------------------------------

    double mn = std::max(lmd.min(), pa.zerotol);

    if (lmd.max() / mn > std::pow(pa.cdet, 1.0 / p))
    {
        arma::mat tmp =  restr2Eigenv(lmd.t(), niini, std::pow(pa.cdet, 1.0 / p), pa.zerotol);

        return tmp.t();
    }

    return lmd;
}

Rcpp::List restrSigmaGPCM(arma::cube SigmaB, arma::vec niini, GPCMPars &pa,      //Rcpp::List paList,
                          bool trace=false, 
                          Rcpp::Nullable<arma::vec> lmd_=R_NilValue,
                          Rcpp::Nullable<arma::cube> OMG_=R_NilValue)
{
    //------------------------------------------------------------
    // Dimensions
    //------------------------------------------------------------

    const unsigned int p = SigmaB.n_rows;
    const unsigned int k = SigmaB.n_slices;

//    paList["p"] = (int)p;
//    paList["k"] = (int)k;
//    GPCMPars pa = parsePars(paList);
    pa.p = (int)p;
    pa.k = (int)k;
    
    const char vol   = pa.pars[0];
    const char shape = pa.pars[1];
    const char rot   = pa.pars[2];

    arma::cube Sigma = SigmaB;

    //------------------------------------------------------------
    // Optional arguments
    //------------------------------------------------------------

    arma::vec lmd;
    arma::cube OMG;

    bool haveLmd = lmd_.isNotNull();
    bool haveOMG = OMG_.isNotNull();

    if (haveLmd)
        lmd = Rcpp::as<arma::vec>(lmd_);

    if (haveOMG)
        OMG = Rcpp::as<arma::cube>(OMG_);

    // DEBUG
    //Rcout << "Entering restrSigmaGPCM()" << std::endl;
    //Rcout << "shw" <<shw << "shb" <<shb << std::endl;
    //SigmaB.print("SigmaB");
    //niini.print("niini");

    //------------------------------------------------------------
    // Number of iterations
    //------------------------------------------------------------

    int maxiterDSR = 1;

    if (pa.pars == "EVE" ||
        pa.pars == "VEE" ||
        pa.pars == "VVE" ||
        pa.pars == "VVV" ||
        pa.pars == "VEV" ||
        pa.pars == "VVI" ||
        pa.pars == "VEI")
    {
        maxiterDSR = pa.maxiterDSR;
    }

    // Parameters not set by the user
    pa.sortsh = 0;
    if(rot == 'E' || rot == 'I') {
        pa.sortsh = 1;
    }

    if(vol == 'E') 
        pa.cdet = 1;
    
    // If OMG is identity, shape restriction parameter within groups is set to 1
    if(shape == 'I') 
        pa.shw = 1;
    
    // if Equal shape is imposed shape restriction parameter between groups is set to 1
    if(shape == 'E') 
        pa.shb = 1;

    // DEBUG
    if(trace)
        Rprintf("\nPARAMETERS: maxiterDSR: %d, maxiterR: %d, maxiterS: %d, cdet: %f, shb: %f, shw: %f\n", maxiterDSR, pa.maxiterR, pa.maxiterS, pa.cdet, pa.shb, pa.shw);

    //------------------------------------------------------------
    // Initialization
    //------------------------------------------------------------

    arma::cube Wk;
    arma::vec wk;

    if (rot == 'E')
    {
        if (!haveLmd || !haveOMG)
        {
            Rcpp::List tmp = initR(SigmaB, niini, p, k, pa.pars);

            lmd = Rcpp::as<arma::vec>(tmp["lmd"]);
            OMG = Rcpp::as<arma::cube>(tmp["Omega"]);
        }

        if (pa.pars == "VVE" || pa.pars == "EVE")
        {
            Wk.set_size(p,p,k);
            wk.set_size(k);

            double sumni = arma::sum(niini);

            for(unsigned int j=0;j<k;++j)
            {
                Wk.slice(j) =
                    (niini(j)/sumni) *
                    SigmaB.slice(j);

                arma::vec eigval;

                arma::eig_sym(eigval,Wk.slice(j));

                wk(j)=eigval.max();
            }
        }
    }
    else if(rot=='V')
    {
        if(!haveLmd || !haveOMG)
        {
            lmd.set_size(k);
            lmd.ones(k);
            OMG.set_size(p,p,k);

            for(unsigned int j=0;j<k;++j)
            {
                arma::vec eigval;
                arma::mat eigvec;

                arma::eig_sym(eigval,eigvec,SigmaB.slice(j));

                if(vol == 'V') {
                    lmd(j)=std::pow(arma::prod(eigval),1.0/p);
                }

                OMG.slice(j)=arma::fliplr(eigvec);
            }
        }
    }
    else
    {
        if(!haveLmd || !haveOMG)
        {
            lmd.ones(k);

            OMG.set_size(p,p,k);

            for(unsigned int j=0;j<k;++j)
            {
                OMG.slice(j).eye();

                if(vol == 'V') {
                    lmd(j)=std::pow(arma::det(SigmaB.slice(j)),1.0/p);
                    if(lmd(j) == 0)
                        lmd(j) = 1;
                }
            }
        }
    }

//  DEBUG
     //Rcout << "Initialization ready" << std::endl;
     //lmd.print("lmd");
     //OMG.print("OMG");

    //------------------------------------------------------------
    // Initial values
    //------------------------------------------------------------

    arma::mat GAM(p,k,arma::fill::ones);

    arma::mat OMGold = OMG.slice(0);

    arma::vec GAMold(p*k);
    arma::mat GAMfc(p,k,arma::fill::zeros);
    GAMold.fill(9999.0);

    arma::vec lmdold(k);
    lmdold.fill(999.0);

    //------------------------------------------------------------
    // Main loop
    //------------------------------------------------------------

    double diffglob = std::numeric_limits<double>::infinity();

    int iter = 0;

    // DEBUG
    // cat("Iter diff_lmd diff_GAM diff_OMG\n")
    if(trace)
        Rcout << "Iter diff_lmd diff_GAM diff_OMG" << std::endl;

    while(diffglob > pa.tolDSR && iter < maxiterDSR)
    {
        ++iter;

        double diffOMG = 0.0;

        //--------------------------------------------------------
        // Update rotation
        //--------------------------------------------------------

        if(iter>1 && rot=='E') {
        
            if(pa.pars=="VVE" || pa.pars=="EVE") {
                Rcpp::List tmp = cpcV(lmd, GAM, OMG.slice(0), Wk, wk, pa.tolR, pa.maxiterR);
                OMG = Rcpp::as<arma::cube>(tmp["Omega"]);
            } else if(pa.pars != "EEE") {
                Rcpp::List tmp = cpcE(lmd, SigmaB, niini);
                OMG = Rcpp::as<arma::cube>(tmp["Omega"]);
            }

            arma::mat OMGnew = OMG.slice(0);
            arma::mat A = OMGnew.t()*OMGold;

            diffOMG = std::abs((double(p)-arma::accu(arma::square(A)))/double(p));
            OMGold = OMGnew;
        }

        //--------------------------------------------------------
        // Update shape
        //--------------------------------------------------------

        GAM = restrshapeGPCM(lmd, OMG, SigmaB, niini, pa);

        // DEBUG
        //  Rcout << "Update GAM (after restrshapeGPCM before sorting):" << std::endl;
        //  GAM.print("GAM");
        //  Rcout << std::endl;

        if(pa.sortsh)
        {
            for(unsigned int j=0;j<k;++j)
                GAMfc.col(j)=arma::sort(GAM.col(j),"descend");
        }else {
            GAMfc = GAM;
        }

       
        //  Rcout << "Update GAM (after restrshapeGPCM after sorting ):" << std::endl;
        //  GAM.print("GAM");
        //  Rcout << std::endl;

        arma::vec GAMnew = arma::vectorise(GAMfc);

        double diffGAM =
            arma::accu(arma::square(GAMnew-GAMold))
            /
            arma::accu(arma::square(GAMold));

        GAMold = GAMnew;

        //--------------------------------------------------------
        // Update determinants
        //--------------------------------------------------------

        lmd =
            restrdeterGPCM(GAM, OMG, SigmaB, niini, pa);

        double difflmd =
            arma::accu(arma::square(lmd-lmdold))
            /
            arma::accu(arma::square(lmdold));

        lmdold = lmd;

        //--------------------------------------------------------
        // Global convergence
        //--------------------------------------------------------

        diffglob =
            std::max(std::max(diffOMG,diffGAM), difflmd);
            
        // cat(iter, difflmd, diffGAM, diffOMG, "\n")
        // Rcout << iter << "   " << difflmd << "   " << diffGAM << "   " << diffOMG << std::endl;
        
        if(trace)
            Rprintf("%d  %f  %f  %f\n", iter, difflmd, diffGAM, diffOMG);            
    }

    //------------------------------------------------------------
    // Reconstruct covariance matrices
    //------------------------------------------------------------

    bool ok = (!GAM.has_nan()) && (GAM.max() > pa.zerotol);

    if(ok)
        Sigma = reconstructSigma(lmd,OMG,GAM);

    return Rcpp::List::create(
        Rcpp::Named("Sigma") = Sigma,
        Rcpp::Named("lmd")   = lmd,
        Rcpp::Named("OMG")   = OMG,
        Rcpp::Named("GAM")   = GAM);
}

// [[Rcpp::export]]
Rcpp::List tclust_restrSigmaGPCM(arma::cube SigmaB, arma::vec niini, Rcpp::List paList, bool trace=false) {
    GPCMPars pa = parsePars(paList);
    return restrSigmaGPCM(SigmaB, niini, pa, trace);
}

