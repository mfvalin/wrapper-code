MOD = gnu
EXT =
OPT = -O2
SHELL = /bin/bash
VAR = TT
DEF =

.SUFFIXES :
.SUFFIXES : .F90 .c .o

OBJECTS = lorenzo_f.o lorenzo_c.o cdf97.o cdf53.o idwt53.o analyze_8x8.o \
	  average_by_2_4_8.o average_2x2.o smooth124.o bi_endian_pack.o \
          misc_timers.o misc_pack.o float_info.o misc_swizzle.o misc_endian.o \
          ieee_quantize.o misc_floats.o verify.o misc_operators.o

.F90.o:
	(module load $(MOD)$(EXT) >/dev/null ; mkdir -p obj_$(MOD) include_$(MOD) ; \
	$$FTN -c -I. -Iinclude_$(MOD) $(OPT) $(DEF) $< -o obj_$(MOD)/$@ ; \
	mv *.mod include_$(MOD)/. )
.c.o:
	(module load $(MOD)$(EXT) >/dev/null ; mkdir -p obj_$(MOD) ; \
	$$CC -c -I. $(OPT) $(DEF) $< -o obj_$(MOD)/$@ ; )

libmisc_$(MOD).a: $(OBJECTS)
	ar rcv $@ obj_$(MOD)/*.o
	mkdir -p include_$(MOD)

all_tests : test_compress test_lorenzo test_pack test_pack_macros \
            test_quantize test_float_info test_helpers test_averages \
            test_cdf97 test_cdf53 test_endian test_endian_f \
            test_ieeequantize

test_compress : test_compress.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; rm -f sample_compressed ; ./$< sample $(VAR) sample_compressed) || true
test_compress.$(MOD).exe : tests/test_compress.F90 libmisc_$(MOD).a
	(module load $(MOD)$(EXT) >/dev/null ; echo "LD_RUN_PATH=$$LD_RUN_PATH" ; \
	$$FTN -fopenmp -I. -Iinclude_$(MOD) $(OPT) $(DEF) $< -L. -lmisc_$(MOD) $$LIBECCC/lib/librmn.a -o $@ -Wl,-rpath,$$LD_RUN_PATH; )

test_lorenzo : test_lorenzo.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_lorenzo.$(MOD).exe : tests/test_lorenzo.F90 libmisc_$(MOD).a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$FTN -I. -Iinclude_$(MOD) $(OPT) $(DEF) $< -L. -lmisc_$(MOD) -o $@ ; )

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
test_bi_endian_pack.$(MOD).exe : tests/test_bi_endian_pack.c libmisc_$(MOD).a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -L. -lmisc_$(MOD) -o $@ )

test_quantize : test_quantize.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_quantize.$(MOD).exe : tests/test_quantize.c libmisc_$(MOD).a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -L. -lmisc_$(MOD) -o $@ )

test_ieeequantize : test_ieeequantize.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_ieeequantize.$(MOD).exe : tests/test_ieeequantize.c libmisc_$(MOD).a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -L. -lmisc_$(MOD) -o $@ )

test_float_info : test_float_info.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_float_info.$(MOD).exe : tests/test_float_info.c
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -L. -lmisc_$(MOD) -o $@ )

test_helpers : test_helpers.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_helpers.$(MOD).exe : tests/test_helpers.c
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -o $@ )

test_average_2x2 : test_average_2x2.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< 128 128 ) || true
test_average_2x2.$(MOD).exe : tests/test_average_2x2.c libmisc_$(MOD).a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -L. -lmisc_$(MOD) -o $@ )

test_averages : test_averages.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_averages.$(MOD).exe : tests/test_averages.c libmisc_$(MOD).a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -L. -lmisc_$(MOD) -o $@ )

test_cdf97 : test_cdf97.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_cdf97.$(MOD).exe : tests/test_cdf.c  libmisc_$(MOD).a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -DCDF97 -I. $^ -L. -lmisc_$(MOD) -o $@ )

test_cdf53 : test_cdf53.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_cdf53.$(MOD).exe : tests/test_cdf.c  libmisc_$(MOD).a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -DCDF53 -I. $^ -L. -lmisc_$(MOD) -o $@ )

test_idwt53 : test_idwt53.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_idwt53.$(MOD).exe : tests/test_idwt.c  libmisc_$(MOD).a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -L. -lmisc_$(MOD) -lm -o $@ )

test_endian : test_endian.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_endian.$(MOD).exe : tests/test_endian.c libmisc_$(MOD).a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -DCDF53 -I. $^ -L. -lmisc_$(MOD) -o $@ )

test_endian_f : test_endian_f.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_endian_f.$(MOD).exe : tests/test_endian_f.F90 libmisc_$(MOD).a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$FTN $(OPT) $(DEF) -DCDF53 -I. $^ -L. -lmisc_$(MOD) -o $@ )

test_floats : test_floats.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_floats.$(MOD).exe : tests/test_floats.c libmisc_$(MOD).a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -DCDF53 -I. $^ -L. -lmisc_$(MOD) -o $@ )

test_smooth124: test_smooth124.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_smooth124.$(MOD).exe : tests/test_smooth124.c libmisc_$(MOD).a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -L. -lmisc_$(MOD) -o $@ )

test_lorenzo_c: test_lorenzo_c.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_lorenzo_c.$(MOD).exe : tests/test_lorenzo_c.c libmisc_$(MOD).a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$CC $(OPT) $(DEF) -I. $^ -L. -lmisc_$(MOD) -o $@ )

clean:
	rm -f *.o *.mod a*.out *.Abs *.exe *.s *.ll *.gch
veryclean: clean
	rm -f libmisc*.a include*/*.mod obj*/*.o

