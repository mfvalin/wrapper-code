MOD = gnu
EXT =
OPT = -O2
SHELL = /bin/bash
VAR = TT
DEF =

.SUFFIXES :
.SUFFIXES : .F90 .c .o

OBJECTS = lorenzo_f.o lorenzo_c.o cdf97.o cdf53.o average_by_2_4_8.o analyze_8x8.o \
          misc_timers.o quantize.o float_info.o misc_swizzle.o misc_endian.o

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
            test_cdf97 test_cdf53 test_endian test_endian_f

test_compress : test_compress.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; rm -f sample_compressed ; ./$< sample $(VAR) sample_compressed) || true
test_compress.$(MOD).exe : tests/test_compress.F90 libmisc_$(MOD).a
	(module load $(MOD)$(EXT) >/dev/null ; \
	$$FTN -I. -Iinclude_$(MOD) $(OPT) $(DEF) $< -L. -lmisc_$(MOD) $$LIBECCC/lib/librmn.a -o $@ ; )

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

test_quantize : test_quantize.$(MOD).exe
	[[ "${EXEC}" != "" ]] && (module load $(MOD)$(EXT) >/dev/null ; ./$< ) || true
test_quantize.$(MOD).exe : tests/test_quantize.c libmisc_$(MOD).a
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

clean:
	rm -f *.o *.mod a*.out *.Abs *.exe *.s
veryclean: clean
	rm -f libmisc*.a include*/*.mod obj*/*.o

