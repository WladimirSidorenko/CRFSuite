# $Id$
TEST_EXTENSIONS = .test
TEST_LOG_DRIVER = env TOP_BUILD_PREFIX='$(top_build_prefix)' \
	TOP_SRCDIR='$(top_srcdir)' AM_TAP_AWK='$(AWK)' \
	$(SHELL) $(top_srcdir)/tap-driver.sh
TESTS = test_sm_1.test

EXTRA_DIST = $(TESTS) \
	test_sm_1_1.input \
        test_sm_1_1.expected

mostlyclean-local: mostlyclean-local-check

.PHONY: mostlyclean-local-check

mostlyclean-local-check:
	-rm -f *.model *.output
