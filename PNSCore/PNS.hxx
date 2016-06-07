#ifndef __PNS_TXX__
#define __PNS_TXX__

#include "PNS.h"

namespace statismo {
    template <typename T> 
        MatrixXd PNS<T>::computeRotMat( const VectorXd& vec ) const {

            MatrixXd rotMat( vec.size(), vec.size() );
            rotMat.setIdentity();

            VectorXd northPole = VectorXd::Zero( vec.size() );
            northPole( vec.size() - 1 ) = 1;

            VectorXd vecUnit = vec;
            vecUnit.normalize();

            double mult = northPole.adjoint() * vecUnit;
            double alpha = acos( mult );

            if ( fabs( mult - 1) < 1e-15 ) {
                return rotMat;
            }
            if ( fabs( mult +1) < 1e-15 ) {
                rotMat*= -1;
                return rotMat;
            }
            VectorXd auxVec = vecUnit - mult*northPole;
            auxVec.normalize();
            MatrixXd auxMat = northPole*auxVec.transpose() - auxVec*northPole.transpose();
            rotMat += sin(alpha)*auxMat + (cos(alpha) - 1)*(northPole*northPole.transpose() + auxVec*auxVec.transpose());

            return rotMat;
        }

    template <typename T>
        MatrixXd PNS<T>::computeRiemannianExpMap( const MatrixXd& mat ) const {

            MatrixXd result( mat.rows()+1, mat.cols() );
            RowVectorXd normVec     = mat.colwise().norm();
            RowVectorXd normVec2    = (normVec.array() < 1e-16).select(0,normVec);
            RowVectorXd sineVec     = normVec2.array().sin();
            RowVectorXd cosVec      = normVec2.array().cos();
            RowVectorXd denoVec     = (normVec2.array() < 1e-16).select(1, normVec);

            RowVectorXd auxVec      = sineVec.array() / (denoVec.array());
            MatrixXd    auxMat      = auxVec.replicate( mat.rows(), 1 );
            auxMat                  = (auxMat.array() * mat.array()).matrix();
            result << auxMat, cosVec;
            return result;
        }
    template <typename T>
        MatrixXd PNS<T>::computeRiemannianLogMap( const MatrixXd& mat ) const {

            MatrixXd result(mat.rows()-1, mat.cols());

            RowVectorXd lastRow = mat.row(mat.rows() - 1);
            //RowVectorXd auxV1   = (lastRow.array() * lastRow.array()).matrix();
            RowVectorXd auxV1   = ( lastRow.array().square() ).matrix();
            RowVectorXd denomV  = ( (1 - auxV1.array() ).sqrt()).matrix();
            RowVectorXd numerV  = ( lastRow.array().acos() ).matrix();
            RowVectorXd auxV3   = ( numerV.array() ) / ( denomV.array() );
            // This line is to check for NaN occurrences.
            // If any entry in auxV3 is NaN then replace it with 1.
            RowVectorXd scale   = (denomV.array() < 1e-64 && numerV.array() < 1e-64).select(1,auxV3); 
            MatrixXd    auxM1   = scale.replicate( mat.rows() - 1, 1 );
            MatrixXd    auxM2   = (mat.topRows( mat.rows() - 1 ));
            result  = (auxM1.array() * auxM2.array()).matrix();

            return result;
        }
    template <typename T>
        double PNS<T>::computeGeodesicMeanS1( const VectorXd& angles ) {
            VectorXd meanCandidate( angles.size() );
            VectorXd auxV1( angles.size() );
            VectorXd auxV2( angles.size() );
            VectorXd theta( angles.size() );
            MatrixXd distMatrix( angles.size(), 3 );
            VectorXd geodVariance( angles.size() );
            double currCandidate(0);
            double currGeodVariance(0);
            int idxToGeodMean(0);
            double geodMean(0);


            // same as theta = mod( angles, 2*pi ) in MATLAB
            theta = angles.unaryExpr( std::ptr_fun(modBy2PI) );
            // Generating mean candidates
            auxV1.setLinSpaced(angles.size(), 0, angles.size()-1);
            auxV1 = (auxV1.array() / auxV1.size()).matrix();

            //meanCandidate = (angles.mean() + 2*PI*auxV1.array()).matrix();
            auxV2 = (angles.mean() + 2*PI*auxV1.array()).matrix();
            //std::cout << "Mean candidate before modding:\n" << auxV2 << std::endl;
            meanCandidate = auxV2.unaryExpr( std::ptr_fun(modBy2PI) );
            //meanCandidate.unaryExpr( std::ptr_fun(modBy2PI));
            //std::cout << "Mean candidate after modding:\n" << meanCandidate << std::endl;
            for(int i=0; i<angles.size(); ++i) {
                double currCandidate = meanCandidate(i);
                distMatrix.col(0) = ((theta.array() - currCandidate).square()).matrix();
                distMatrix.col(1) = ((theta.array() - currCandidate + 2*PI).square()).matrix();
                distMatrix.col(2) = ((currCandidate - theta.array() + 2*PI).square()).matrix();
                currGeodVariance = (distMatrix.rowwise().minCoeff()).sum();
                geodVariance(i) = currGeodVariance;
            }
            double minGeodVarianceValue = geodVariance.minCoeff(&idxToGeodMean);
            geodMean = modBy2PI( meanCandidate( idxToGeodMean ) );
        }
    template <typename T>
        double PNS<T>::LMsphereFit( const MatrixXd& data, const itype flag, VectorXd& x ) const {
            // returns radius and x is the result of the optimization
            int info;
            int sizeX = data.rows();
            int sizeY = data.cols();
            double r = -1;
            assert( data.rows() == x.size() ); 
            // create functor
            sphereResFunctor functor(sizeX,sizeY, data, flag);
            LevenbergMarquardt<sphereResFunctor> optimizer(functor);
            info = optimizer.minimize(x);
            if ( info > 0 ) {
                VectorXd fvec( data.cols() );
                r = functor(x, fvec);
            }
            return r;
        }
    template <typename T>
        double PNS<T>::objFn( const VectorXd& center, const MatrixXd& data, const double r ) const {
            double result = -1;
            // take a dot product between center and the data
            VectorXd aux1 = center.adjoint()*data;
            VectorXd aux2 = ( ( aux1.array().acos() - r ).square() ).matrix();
            result = aux2.array().mean();
            return result;
        }

    template <typename T> double
        PNS<T>::computeSubSphere( const MatrixXd& data, const itype flag, VectorXd& center, double& error ) const {
            // set up initial parameters
            int trial = 0;
            int nTrial = 100;
            error = 1e+10;
            double curr_error = 0;
            double error_diff = fabs( curr_error - error );
            double TOL = 1e-10;
            double radius = -1;

            MatrixXd rotMat( center.size(), center.size() );
            MatrixXd rotatedData( data.rows(), data.cols() );
            MatrixXd TanProjData( data.rows()-1, data.cols() );
            VectorXd tempCenter( center.size() -1 );

            VectorXd newCenter( center.size() );

            while ( (error_diff > TOL) && (trial < nTrial) ) {
                center.normalize();
                rotMat = computeRotMat( center );
                rotatedData = rotMat*data;

                TanProjData = computeRiemannianLogMap( rotatedData );
                tempCenter.setZero();
                radius = LMsphereFit( TanProjData, flag, tempCenter );
                assert( radius > 0 ); // for sanity check

                if ( radius > PI ) {
                    radius = PI/2;
                    JacobiSVD<MatrixXd> svd( TanProjData, Eigen::ComputeThinU );
                    tempCenter = svd.matrixU().rightCols(1).array() * (PI/2);
                }
                newCenter = computeRiemannianExpMap( tempCenter );
                center = rotMat.inverse() * newCenter;

                curr_error = objFn( center, data, radius );
                error_diff = fabs( curr_error - error );
                error = curr_error;
                trial++;
            }

            return radius;
        }
    template <typename T>
        double PNS<T>::getSubSphere( const MatrixXd& data, const itype flag, VectorXd& center) const {
            // first svd
            // There is the sign ambiguity in svd module
            // I am not sure if this will affect the end result
            // Inputs for svd
            MatrixXd X = data;  // Input data is in column major matrix
            X.transposeInPlace(); // So transpose it to be used for svd and pca
            double r1 = -1;
            double err1 = 1e+10;

            JacobiSVD<MatrixXd> svd( data, ComputeThinU );
            VectorXd center1( svd.matrixU().rightCols(1) );
            r1 = computeSubSphere(data, flag, center1, err1);

            // Inputs for pca
            double r2 = -1;
            double err2 = 1e+10;

            unsigned n = X.rows();

            RowVectorXd mu = X.colwise().mean(); // needs to be row vector
            MatrixXd X0 = X.rowwise() - mu; // center data
            JacobiSVD<MatrixXd> pca( X0.transpose()*X0*(1.0/(n-1) ), ComputeThinU );
            VectorXd center2( pca.matrixU().rightCols(1) );
            r2 = computeSubSphere(data, flag, center2, err2);
            if ( err1 < err2 ) {
                center = center1;
                return r1;
            }
            else {
                center = center2;
                return r2;
            }
        }
    template <typename T>
        MatrixXd PNS<T>::compute( const double& alpha, const double& R ) {
            // returns residual matrix of size (n-1) x p if the data_ is n x p
            // also saves PNS data structures internally
            // This includes the following:
            // 1. mean
            // 2. radii of each nested sphere
            // 3. orthoaxis
            // 4. distance
            // 5. p-values (TO be implemented later)


            unsigned int k = data_.rows();
            unsigned int n = data_.cols();

            JacobiSVD<MatrixXd> svd( data_, ComputeThinU );
            unsigned int maxd = (svd.singularValues().array() < 1e-16).count();
            // there is some maxd stuff that I will do later
            if ( (maxd == 0) || (k > n) ) {
                maxd = std::min(k,n) + 1;
            }
            unsigned int nullspdim = k - maxd + 1; // dimension of subspace that contains no data
            unsigned int dm = maxd - 2; // intrinsic dimension of the smallested nested sphere that contains variation, i.e., the dimension to whcih can be trivailly reduced
            if ( k == 3 ) {
                nullspdim = 0;
                dm = 2;
            }
            MatrixXd currentSubSphere;
            if (nullspdim > 0) {
                currentSubSphere = (svd.matrixU().leftCols(dm+1)).adjoint()*data_;
            }
            else {
                currentSubSphere = data_;
            }
            // allocating residual matrix
            MatrixXd residualMatrix( dm , data_.cols() );
            // TODO: check to see if I can get away without creating intermediate data
            MatrixXd nestedSphere;
            MatrixXd currRotation;
            MatrixXd auxM1;
            MatrixXd auxM2;


            VectorXd center;
            VectorXd currResidual;
            RowVectorXd auxV1;
            double r = -1;
            switch (flag_) {
                //case SEQTEST:
                // Need to call getSubSphere twice with one with Great Circle option and Small Circle Option
                //    break;
                case SMALL: case GREAT:
                    for(int i = 0; i < dm-1; ++i ) {
                        // First Estimate the best fitting subsphere
                        // Either with SMALL sphere if itype is SMALL
                        // or with GREAT sphere if itype is GREAT
                        r = getSubSphere(currentSubSphere, flag_, center);
                        currResidual = ((center.transpose()*currentSubSphere).array()).acos() - r;
                        // save subsphere parameters
                        orthaxis_.push_back(center);
                        dist_.push_back(r);
                        // record residual
                        residualMatrix.row(i) = currResidual;
                        // projection to subsphere and transformation to isomorphic sphere
                        currRotation = computeRotMat(center);
                        // To make sure nestedSphere has the same size as currentSubSphere
                        nestedSphere.resize( currentSubSphere.rows(), currentSubSphere.cols() );
                        nestedSphere = currRotation*currentSubSphere;

                        auxM1 = nestedSphere.topRows(dm-i);
                        auxV1 = ( 1 - nestedSphere.row( dm-i ).array().square() ).sqrt();
                        auxM2 = auxV1.replicate( dm - i, 1 );
                        // reduce the size of currentSubSphere
                        currentSubSphere.resize(auxM1.rows(), auxM1.cols());
                        currentSubSphere = auxM1.array() / auxM2.array();
                    }
                    break;
                default:
                    break;
            }
            // currentSubSphere has intrinsic dimension 1
            
            VectorXd S1toRadian = currentSubSphere.row(1).binaryExpr( currentSubSphere.row(0), MakeAtan2Op<double>() ); 
            double meantheta = -1;
            meantheta = computeGeodesicMeanS1( S1toRadian );
            VectorXd meantheta_vector(1);
            meantheta_vector(0) = meantheta;
            orthaxis_.push_back( meantheta_vector );
            VectorXd auxV2 = S1toRadian.array() - meantheta + PI;
            VectorXd auxV3 = auxV2.unaryExpr( std::ptr_fun(modBy2PI) );
            residualMatrix.bottomRows(1) = (auxV3.array() - PI).transpose();
            // computing radii
            Eigen::Map< VectorXd > aux_map_vector_dist( dist_.data(), dist_.size() );
            double result = 1;
            radii_.push_back(result);
            for(int i = 0; i < dm-1; ++i) {
                result = aux_map_vector_dist.head(i+1).array().sin().prod();
                radii_.push_back( result );
            }

            if (nullspdim > 0) {
                basisu_flag_ = true;
                basisu_ = svd.matrixU().leftCols(dm+1);
            }
            Eigen::Map < VectorXd > aux_map_vector_radii( radii_.data(), radii_.size() );
            residualMatrix = (( aux_map_vector_radii.replicate( 1, n ).array() * residualMatrix.array() ).colwise().reverse()).eval();
            return residualMatrix;

        }
    template <typename T>
        MatrixXd PNS<T>::computeS2E( const MatrixXd& sphereData ) const {
            //TODO: add a dimensionality check
            MatrixXd result;
            MatrixXd current_sphere;
            if(basisu_flag_) {
                current_sphere = basisu_*sphereData;
            }
            else {
                current_sphere = sphereData;
            }

            unsigned int n_dim = current_sphere.rows();
            unsigned int n_case = current_sphere.cols();

            MatrixXd res_matrix = MatrixXd::Zero( n_dim-1, n_case );
            MatrixXd nested_sphere;
            MatrixXd aux_scale_matrix;

            VectorXd current_orthaxis;
            
            RowVectorXd current_residual(n_case);

            double current_radius;

            for(int i=0; i<n_dim-2; ++i) { 
                current_orthaxis = orthaxis_[i];
                current_radius = dist_[i];
                current_residual = (current_orthaxis.transpose()*current_sphere).array().acos() - current_radius;
                res_matrix.row(i) = current_residual;
                nested_sphere = computeRotMat(current_orthaxis)*current_sphere;
                aux_scale_matrix = ( 1 - nested_sphere.bottomRows(1).array().square() ).sqrt().replicate(n_dim-i-1, 1 );
                current_sphere = nested_sphere.topRows(n_dim-i-1).array() / aux_scale_matrix.array();
            }
            
            VectorXd S1toRadian = current_sphere.row(1).binaryExpr( current_sphere.row(0), MakeAtan2Op<double>() ); 
            double geomean = orthaxis_.back()[0];
            //theta = angles.unaryExpr( std::ptr_fun(modBy2PI) );
            VectorXd devS1 = ( S1toRadian.array() - geomean + PI ).unaryExpr( std::ptr_fun(modBy2PI) ) - PI;
            res_matrix.row( n_dim - 2 ) = devS1;
            Eigen::Map< const VectorXd > aux_mapped_radii_vector( radii_.data(), radii_.size() );
            result = (aux_mapped_radii_vector.replicate(1,n_case).array() * res_matrix.array() ).colwise().reverse();
            return result;
        }
    template <typename T>
        MatrixXd PNS<T>::computeE2S( const MatrixXd& euclideanData ) const {
            MatrixXd result;
            unsigned int n_dim = euclideanData.rows();
            unsigned int n_case = euclideanData.cols();

            result.resize(2,n_case);

            vector< VectorXd > nestedsphere_orthaxis( orthaxis_.rbegin()+1, orthaxis_.rend());
            vector< double > nestedsphere_radius( dist_.rbegin(), dist_.rend() );
            //VectorXd geomean = orthaxis_.back();
            double geomean = orthaxis_.back()[0];
            // standardize the coordinates
            vector< double > nestedsphere_radii( radii_.rbegin(), radii_.rend() );
            Eigen::Map< VectorXd >  aux_radii_vector( nestedsphere_radii.data(), nestedsphere_radii.size() );
            MatrixXd res_standardized = euclideanData.array() / ( aux_radii_vector.replicate( 1, n_case ) ).array();
            // to S^1
            if ( n_dim > 0 ) {
                result.row(0) = (res_standardized.row(0).array() + geomean).cos();
                result.row(1) = (res_standardized.row(0).array() + geomean).sin();
            }
            // S^1 to S^2
            if ( n_dim > 1 ) {
                cout << "Orthoaxis:\n" << nestedsphere_orthaxis[0] << endl;
                MatrixXd aux_rotation_matrix = computeRotMat( nestedsphere_orthaxis[0] );
                aux_rotation_matrix.transposeInPlace();

                MatrixXd aux_matrix_scale = ( (res_standardized.row(1).array() + nestedsphere_radius[0] ).sin()).replicate(2,1);
                result = aux_matrix_scale.array()*result.array();

                result.conservativeResize(3, n_case);
                result.row(2) = (res_standardized.row(1).array() + nestedsphere_radius[0]).cos();

                result = aux_rotation_matrix*result;
            }
            // S^2 to S^d
            if ( n_dim > 2 ) {
                // be careful with the dimension and number of iteration
                for(int i=0;i<n_dim-2;++i) {
                    MatrixXd aux_rotation_matrix_2 = computeRotMat( nestedsphere_orthaxis[i+1] );
                    aux_rotation_matrix_2.transposeInPlace();

                    MatrixXd aux_matrix_scale_2 = ( (res_standardized.row(i+2).array() + nestedsphere_radius[i+1] ).sin()).replicate(result.rows(),1);
                    result = aux_matrix_scale_2.array()*result.array();

                    result.conservativeResize(result.rows()+1, n_case);
                    result.bottomRows(1) = (res_standardized.row(i+2).array() + nestedsphere_radius[i+1]).cos();

                    result = aux_rotation_matrix_2*result;
                }
            }
            if ( basisu_flag_ ) 
                result = (basisu_*result);
            return result;
        }

} // namespace statismo
#endif // __PNS_TXX__
