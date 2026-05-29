# -*-Makefile-*-

all: McLachlan_ADM.out McLachlan_BSSN.out McLachlan_BSSN_bench.out McLachlan_BSSN_bench8.out McLachlan_BSSN_CL.out McLachlan_BSSN_ND.out McLachlan_BSSN_NV.out McLachlan_CCZ4.out McLachlan_ADMConstraints.out McLachlan_ADMQuantities.out WaveToy.out WaveToyFO.out hydro.out Carlile_BSSN.out
	@echo
	@echo "The Cactus thorns are up to date."
	@echo

McLachlan_ADM.out: McLachlan_ADM.m
	rm -rf ML_ADM
	./runmath.sh $^
	for thorn in ML_ADM; do					\
		./copy-if-changed.sh $$thorn ../$$thorn;	\
	done

McLachlan_BSSN.out: McLachlan_BSSN.m
	rm -rf ML_BSSN ML_BSSN_Helper
	env ML_CODE=BSSN ./runmath.sh $^ McLachlan_BSSN
	./create-helper-thorn.sh ML_BSSN ML_BSSN
	./copy-if-changed.sh ML_BSSN ../ML_BSSN
	./copy-if-changed.sh ML_BSSN_Helper ../ML_BSSN_Helper

McLachlan_BSSN_bench.out: McLachlan_BSSN.m
	rm -rf ML_BSSN_bench ML_BSSN_bench_Helper
	env ML_CODE=BSSN_bench ./runmath.sh $^ McLachlan_BSSN_bench
	./create-helper-thorn.sh ML_BSSN ML_BSSN_bench
	./copy-if-changed.sh ML_BSSN_bench ../ML_BSSN_bench
	./copy-if-changed.sh ML_BSSN_bench_Helper ../ML_BSSN_bench_Helper

McLachlan_BSSN_bench8.out: McLachlan_BSSN.m
	rm -rf ML_BSSN_bench8 ML_BSSN_bench8_Helper
	env ML_CODE=BSSN_bench8 ./runmath.sh $^ McLachlan_BSSN_bench8
	./create-helper-thorn.sh ML_BSSN ML_BSSN_bench8
	./copy-if-changed.sh ML_BSSN_bench8 ../ML_BSSN_bench8
	./copy-if-changed.sh ML_BSSN_bench8_Helper ../ML_BSSN_bench8_Helper

McLachlan_BSSN_CL.out: McLachlan_BSSN.m
	rm -rf ML_BSSN_CL ML_BSSN_CL_Helper
	env ML_CODE=BSSN_CL ./runmath.sh $^ McLachlan_BSSN_CL
	./create-helper-thorn.sh ML_BSSN ML_BSSN_CL
	./copy-if-changed.sh ML_BSSN_CL ../ML_BSSN_CL
	./copy-if-changed.sh ML_BSSN_CL_Helper ../ML_BSSN_CL_Helper

McLachlan_BSSN_NV.out: McLachlan_BSSN.m
	rm -rf ML_BSSN_NV ML_BSSN_NV_Helper
	env ML_CODE=BSSN_NV ./runmath.sh $^ McLachlan_BSSN_NV
	./create-helper-thorn.sh ML_BSSN ML_BSSN_NV
	./copy-if-changed.sh ML_BSSN_NV ../ML_BSSN_NV
	./copy-if-changed.sh ML_BSSN_NV_Helper ../ML_BSSN_NV_Helper

McLachlan_BSSN_ND.out: McLachlan_BSSN.m
	rm -rf ML_BSSN_ND ML_BSSN_ND_Helper
	env ML_CODE=BSSN_ND ./runmath.sh $^ McLachlan_BSSN_ND
	./create-helper-thorn.sh ML_BSSN ML_BSSN_ND
	./copy-if-changed.sh ML_BSSN_ND ../ML_BSSN_ND
	./copy-if-changed.sh ML_BSSN_ND_Helper ../ML_BSSN_ND_Helper

McLachlan_CCZ4.out: McLachlan_BSSN.m
	rm -rf ML_CCZ4 ML_CCZ4_Helper
	env ML_CODE=CCZ4 ./runmath.sh $^ McLachlan_CCZ4
	./create-helper-thorn.sh ML_BSSN ML_CCZ4
	./copy-if-changed.sh ML_CCZ4 ../ML_CCZ4
	./copy-if-changed.sh ML_CCZ4_Helper ../ML_CCZ4_Helper

McLachlan_ADMConstraints.out: McLachlan_ADMConstraints.m
	rm -rf ML_ADMConstraints*
	./runmath.sh $^
	for thorn in ML_ADMConstraints*; do			\
		./copy-if-changed.sh $$thorn ../$$thorn;	\
	done

McLachlan_ADMQuantities.out: McLachlan_ADMQuantities.m
	rm -rf ML_ADMQuantities*
	./runmath.sh $^
	for thorn in ML_ADMQuantities*; do			\
		./copy-if-changed.sh $$thorn ../$$thorn;	\
	done

WaveToy.out: WaveToy.m
	rm -rf ML_WaveToy ML_WaveToy_CL
	./runmath.sh $^
	for thorn in ML_WaveToy ML_WaveToy_CL; do		\
		./copy-if-changed.sh $$thorn ../$$thorn;	\
	done

WaveToyFO.out: WaveToyFO.m
	rm -rf ML_WaveToyFO
	./runmath.sh $^
	for thorn in ML_WaveToyFO; do				\
		./copy-if-changed.sh $$thorn ../$$thorn;	\
	done

hydro.out: hydro.m
	rm -rf ML_hydro
	./runmath.sh $^
	for thorn in ML_hydro; do				\
		./copy-if-changed.sh $$thorn ../$$thorn;	\
	done

Carlile_BSSN.out: Carlile_BSSN.m
	rm -rf CL_BSSN*
	./runmath.sh $^
	for thorn in CL_BSSN*; do					   \
		./create-helper-thorn.sh CL_BSSN $$thorn;		   \
		./copy-if-changed.sh $$thorn ../$$thorn;		   \
		./copy-if-changed.sh $${thorn}_Helper ../$${thorn}_Helper; \
	done

clean:
	rm -rf ML_ADM ML_BSSN* ML_CCZ4* ML_ADMConstraints* ML_ADMQuantities*
	rm -rf CL_BSSN*
	rm -rf ML_WaveToy ML_WaveToy_CL ML_WaveToyFO
	rm -rf ML_WaveToyS2 ML_WaveToyMetricS2 ML_MaxwellS2 ML_BSSNS2
	rm -rf ML_hydro
	rm -f McLachlan_ADM.out McLachlan_ADM.err
	rm -f McLachlan_BSSN.out McLachlan_BSSN.err
	rm -f McLachlan_BSSN_bench.out McLachlan_BSSN_bench.err
	rm -f McLachlan_BSSN_bench8.out McLachlan_BSSN_bench8.err
	rm -f McLachlan_BSSN_CL.out McLachlan_BSSN_CL.err
	rm -f McLachlan_BSSN_ND.out McLachlan_BSSN_ND.err
	rm -f McLachlan_BSSN_NV.out McLachlan_BSSN_NV.err
	rm -f McLachlan_CCZ4.out McLachlan_CCZ4.err
	rm -f McLachlan_ADMConstraints.out McLachlan_ADMConstraints.err
	rm -f McLachlan_ADMQuantities.out McLachlan_ADMQuantities.err
	rm -f WaveToy.out WaveToy.err
	rm -f WaveToyFO.out WaveToyFO.err
	rm -f WaveToyS2.out WaveToyS2.err
	rm -f WaveToyMetricS2.out WaveToyMetricS2.err
	rm -f MaxwellS2.out MaxwellS2.err
	rm -f BSSNS2.out BSSNS2.err
	rm -f hydro.out hydro.err
	rm -f Carlile_BSSN.out Carlile_BSSN.err

.PHONY: all clean
