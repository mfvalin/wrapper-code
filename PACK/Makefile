MOD = gnu
EXT =
OPT = -O2
SHELL = /bin/bash
VAR = TT
DEF =
LIBDIR = lib_${MOD}
VPATH = obj_${MOD}:rmn

.SUFFIXES :
.SUFFIXES : .F90 .c .o

# %.o : %.h
%.c.o : %.c
	(module load $(MOD)$(EXT) >/dev/null ; mkdir -p obj_$(MOD) ; \
	$$CC -c -I. -I$$LIBECCC/include $(OPT) $(DEF) $< -o obj_$(MOD)/$@ ; )
%.F90.o: %.F90
	(module load $(MOD)$(EXT) >/dev/null ; mkdir -p obj_$(MOD) include_$(MOD) ; \
	$$FTN -c -I. -Iinclude_$(MOD) -I$$LIBECCC/include $(OPT) $(DEF) $< -o obj_$(MOD)/$@ ; \
	mv *.mod include_$(MOD)/. )

# the order of the Fortran objects is NOT arbitrary, there are dependencies
OBJECTS = lorenzo_c.c.o cdf97.c.o cdf53.c.o idwt53.c.o analyze_8x8.c.o misc_analyze.c.o \
	  average_by_2_4_8.c.o average_2x2.c.o smooth124.c.o bi_endian_pack.c.o \
          misc_timers.c.o misc_pack.c.o float_info.c.o misc_swizzle.c.o misc_endian.c.o \
          ieee_quantize.c.o misc_floats.c.o verify.c.o misc_operators.c.o misc_zfp.c.o \
	  zfpx_quantize.c.o misc_encode.c.o pack_test_utils.c.o \
          lorenzo_f.F90.o misc_operators_mod.F90.o misc_helpers_mod.F90.o rmn_timers_mod.F90.o \
          rmn_stream_pack_mod.F90.o rmn_zfpx_mod.F90.o analyze_data.F90.o rmn_encode_mod.F90.o

INCLUDES = average_2x2.h    cdf53.h        cpu_type.h        idwt53.h \
           lorenzo.h        misc_floats.h  misc_operators.h  misc_swizzle.h  \
           misc_types.h     pack_stream.h  bi_endian_pack.h  cdf97.h  \
           ieee_quantize.h  misc_endian.h  misc_helpers.h    misc_pack.h \
           misc_timers.h    pack_macros.h  smooth124.h       fast_endian.h \
           misc_zfp.h       misc_analyze.h misc_simd.h       misc_encode.h

${LIBDIR}/libmisc.a: recompile_${MOD} $(OBJECTS)
	mkdir -p ${LIBDIR}
	ar rcv $@ obj_$(MOD)/*.o
	mkdir -p include_$(MOD)

recompile_${MOD} :  $(INCLUDES)
	rm -f obj_$(MOD)/*.o
	touch recompile_${MOD}

all_tests : test_compress test_lorenzo test_pack test_pack_macros \
            test_quantize test_float_info test_helpers test_averages \
            test_cdf97 test_cdf53 test_endian test_endian_f \
            test_ieeequantize test_encode

test_pack_value : test_pack_value.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_pack_value.$(MOD).exe : tests/test_pack_value.c ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC -I. -Iinclude_$(MOD) $(OPT) $(DEF) $< -L${LIBDIR} -lmisc -o $@ ; )

test_timers : test_timers.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_timers.$(MOD).exe : tests/test_timers.c ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC -I. -Iinclude_$(MOD) $(OPT) $(DEF) $< -L${LIBDIR} -lmisc -o $@ ; )

test_operators_c : test_operators_c.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_operators_c.$(MOD).exe : tests/test_operators_c.c ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC -I. -Iinclude_$(MOD) $(OPT) $(DEF) $< -L${LIBDIR} -lmisc -lm -o $@ ; )

test_compress_c : test_compress_c.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_compress_c.$(MOD).exe : tests/test_compress_c.c ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC -I. -Iinclude_$(MOD) $(OPT) $(DEF) $< -L${LIBDIR} -lmisc -o $@ ; )

test_compress : test_compress.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; \
	rm -f sample_compressed ; ./$< sample $(VAR) sample_compressed) || true
test_compress.$(MOD).exe : tests/test_compress.F90 ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; echo "LD_RUN_PATH=$$LD_RUN_PATH" ; \
	$$FTN -fopenmp -I. -Iinclude_$(MOD) $(OPT) $(DEF) $< -L${LIBDIR} \
	-lmisc $$LIBECCC/lib/librmn.a $$LIBECCC/lib/libzfp.a -o $@ -Wl,-rpath,$$LD_RUN_PATH; )

test_bits_needed : test_bits_needed.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_bits_needed.$(MOD).exe : tests/test_bits_needed.c ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC -I. -Iinclude_$(MOD) $(OPT) $(DEF) $< -L${LIBDIR} -lmisc -o $@ ; )

test_encode : test_encode.$(MOD).exe misc_encode.h
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_encode.$(MOD).exe : tests/test_encode.c ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC -I. -Iinclude_$(MOD) $(OPT) $(DEF) $< -L${LIBDIR} -lmisc -o $@ ; )

test_operators : test_operators.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_operators.$(MOD).exe : tests/test_operators.F90 ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$FTN -I. -Iinclude_$(MOD) $(OPT) $(DEF) $< -L${LIBDIR} -lmisc -o $@ ; )

test_analyze : test_analyze.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_analyze.$(MOD).exe : tests/test_analyze.F90 ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$FTN -I. -Iinclude_$(MOD) $(OPT) $(DEF) $< -L${LIBDIR} -lmisc -o $@ ; )

test_lorenzo : test_lorenzo.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_lorenzo.$(MOD).exe : tests/test_lorenzo.F90 ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$FTN -I. -Iinclude_$(MOD) $(OPT) $(DEF) $< -L${LIBDIR} -lmisc -o $@ ; )

test_pack : test_pack.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_pack.$(MOD).exe : tests/test_pack.c
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -o $@ )

test_pack_macros : test_pack_macros.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_pack_macros.$(MOD).exe : tests/test_pack_macros.c
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -o $@ )

test_bi_endian_pack : test_bi_endian_pack.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_bi_endian_pack.$(MOD).exe : tests/test_bi_endian_pack.c ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -L${LIBDIR} -lmisc -o $@ )

test_zfp : test_zfp.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< .001 ) || true
test_zfp.$(MOD).exe : tests/test_zfp.c ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. -I$${LIBECCC}/include $^ -L${LIBDIR} -L$${LIBECCC}/lib -lmisc -lzfp -lm -o $@ )

test_quantize : test_quantize.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_quantize.$(MOD).exe : tests/test_quantize.c ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -L${LIBDIR} -lmisc -o $@ )

test_ieeequantize : test_ieeequantize.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_ieeequantize.$(MOD).exe : tests/test_ieeequantize.c ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -L${LIBDIR} -lmisc -o $@ )

test_float_info : test_float_info.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_float_info.$(MOD).exe : tests/test_float_info.c
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -L${LIBDIR} -lmisc -o $@ )

test_helpers : test_helpers.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_helpers.$(MOD).exe : tests/test_helpers.c
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -o $@ )

test_average_2x2 : test_average_2x2.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< 128 128 ) || true
test_average_2x2.$(MOD).exe : tests/test_average_2x2.c ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -L${LIBDIR} -lmisc -o $@ )

test_averages : test_averages.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_averages.$(MOD).exe : tests/test_averages.c ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -L${LIBDIR} -lmisc -o $@ )

test_cdf97 : test_cdf97.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_cdf97.$(MOD).exe : tests/test_cdf.c  ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -DCDF97 -I. $^ -L${LIBDIR} -lmisc -o $@ )

test_cdf53 : test_cdf53.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_cdf53.$(MOD).exe : tests/test_cdf.c  ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -DCDF53 -I. $^ -L${LIBDIR} -lmisc -o $@ )

test_idwt53 : test_idwt53.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_idwt53.$(MOD).exe : tests/test_idwt.c  ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -L${LIBDIR} -lmisc -lm -o $@ )

test_endian : test_endian.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_endian.$(MOD).exe : tests/test_endian.c ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -DCDF53 -I. $^ -L${LIBDIR} -lmisc -o $@ )

test_endian_f : test_endian_f.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_endian_f.$(MOD).exe : tests/test_endian_f.F90 ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$FTN $(OPT) $(DEF) -DCDF53 -I. $^ -L${LIBDIR} -lmisc -o $@ )

test_floats : test_floats.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_floats.$(MOD).exe : tests/test_floats.c ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -DCDF53 -I. $^ -L${LIBDIR} -lmisc -o $@ )

test_smooth124: test_smooth124.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_smooth124.$(MOD).exe : tests/test_smooth124.c ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -L${LIBDIR} -lmisc -o $@ )

test_lorenzo_c: test_lorenzo_c.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_lorenzo_c.$(MOD).exe : tests/test_lorenzo_c.c ${LIBDIR}/libmisc.a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -L${LIBDIR} -lmisc -o $@ )

clean:
	rm -f *.o *.mod a*.out *.Abs *.exe *.s *.ll *.gch recompile_*
veryclean: clean
	rm -f lib_*/libmisc.a include*/*.mod obj*/*.o recompile_*

