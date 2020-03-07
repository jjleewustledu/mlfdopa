classdef Test_Ki_Criswell < matlab.unittest.TestCase
	%% TEST_KI_CRISWELL 

	%  Usage:  >> results = run(mlfdopa_unittest.Test_Ki_Criswell)
 	%          >> result  = run(mlfdopa_unittest.Test_Ki_Criswell, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 11-Dec-2018 14:44:42 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfdopa/test/+mlfdopa_unittest.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
        pwd0
        pwd1 = fullfile(getenv('HOME'), 'MATLAB-Drive', 'mlfdopa', 'data', '')
        csv = 'mars_10722_fdopa1_333_fwhm10_freesurfer_rois.csv'
 		registry
 		testObj
 	end

	methods (Test)
        function test_ctor(this)
            import mlfdopa.*;
            this.verifyEqual(this.testObj.csv, this.csv);
            this.verifyEqual(length(this.testObj.estimates), 12);
            this.verifyEqual(length(this.testObj.estimates), this.testObj.Nregions);
            e = this.testObj.estimates(1);
            this.verifyEqual(e.Ki, 0.008898576350779, 'RelTol', 1e-9);
            this.verifyEqual(e.CIupper, 0.009830319720228, 'RelTol', 1e-9);
            this.verifyEqual(e.CIlower, 0.007966832981330, 'RelTol', 1e-9);
            this.verifyEqual(e.pValue, 7.924070234490734e-13, 'RelTol', 1e-9);
            this.verifyEqual(e.Rsquared, 0.962432782535622, 'RelTol', 1e-9);
        end
        function test_plots(this)
            this.testObj.plotIntermed;
            this.testObj.plotPatlak;
            this.testObj.plotTacs;
        end
        function test_writetable(this)
            T = this.testObj.writetable;
            disp(T);
            this.verifyTrue(lexist( ...
                fullfile(this.pwd1, [mybasename(this.testObj.csv) this.testObj.suffix '.csv']), 'file'));
        end
        function test_Kocc(this)
            [est,obj] = mlfdopa.Ki_Criswell.Kocc(this.csv);
            this.verifyEqual(est{'Kocc',:}, obj.Ki');
            this.verifyEqual(size(est), [5 12]);
            this.verifyEqual(est{'Kocc',1},   0.00889857635077907, 'RelTol', 1e-9);
            this.verifyEqual(est{'Kocc',end}, 0.00341814892888871, 'RelTol', 1e-9);
        end
	end

 	methods (TestClassSetup)
		function setupKi_Criswell(this)
            import mlfdopa.*;
 			this.pwd0 = pushd(this.pwd1);
            this.testObj_ = Ki_Criswell(this.csv);
 			this.addTeardown(@this.cleanTest);
 		end
	end

 	methods (TestMethodSetup)
		function setupKi_CriswellTest(this)
            this.testObj = this.testObj_;
 		end
	end

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
		function cleanTest(this)
            popd(this.pwd0);
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

