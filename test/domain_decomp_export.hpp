
#ifndef LIB_EXPORT_H
#define LIB_EXPORT_H

#ifdef DOMAIN_DECOMP_STATIC_DEFINE
#  define LIB_EXPORT
#  define DOMAIN_DECOMP_NO_EXPORT
#else
#  ifndef LIB_EXPORT
#    ifdef domain_decomp_EXPORTS
        /* We are building this library */
#      define LIB_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define LIB_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef DOMAIN_DECOMP_NO_EXPORT
#    define DOMAIN_DECOMP_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef DOMAIN_DECOMP_DEPRECATED
#  define DOMAIN_DECOMP_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef DOMAIN_DECOMP_DEPRECATED_EXPORT
#  define DOMAIN_DECOMP_DEPRECATED_EXPORT LIB_EXPORT DOMAIN_DECOMP_DEPRECATED
#endif

#ifndef DOMAIN_DECOMP_DEPRECATED_NO_EXPORT
#  define DOMAIN_DECOMP_DEPRECATED_NO_EXPORT DOMAIN_DECOMP_NO_EXPORT DOMAIN_DECOMP_DEPRECATED
#endif

/* NOLINTNEXTLINE(readability-avoid-unconditional-preprocessor-if) */
#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef DOMAIN_DECOMP_NO_DEPRECATED
#    define DOMAIN_DECOMP_NO_DEPRECATED
#  endif
#endif

#endif /* LIB_EXPORT_H */
