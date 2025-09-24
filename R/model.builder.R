#' Helper function that generates the misspecified
#' @keywords internal
#' @param seed see gen_fit2
#' @param sm ...
#' @param mm ...
#' @param rcv ...
#' @param pop.mod ...
#' @param magni see sm.misspec, mm.misspec
#' @param tc see sm.misspec
#' @param te see rcv.misspec
#' @noRd
model.builder <-
  function(seed = NULL,
           sm,
           mm,
           rcv,
           pop.mod,
           magni = NULL,
           tc = NULL,
           te = NULL) {
    #Remove unintended whitespaces first for compatibility
    pop.mod <- stringr::str_replace_all(pop.mod, stringr::fixed(" "), "")
    #000
    if (sm == 0 & mm == 0 & rcv == 0) {
      pop.mod.new <- pop.mod
    }
    #x00
    if (sm > 0 & mm == 0 & rcv == 0) {
      pop.mod.new <- sm.misspec(pop.mod,
                                mc = sm,
                                tc = tc,
                                seed = seed)
    }
    #0x0
    if (sm == 0 & mm > 0 & rcv == 0) {
      pop.mod.new <-
        mm.misspec(pop.mod,
                   ml = mm,
                   magni = magni,
                   seed = seed)
    }
    #00x
    if (sm == 0 & mm == 0 & rcv > 0) {
      pop.mod.new <- rv.misspec(pop.mod,
                                me = rcv,
                                te = te,
                                seed = seed)
    }
    #0xx
    if (sm == 0 & mm > 0 & rcv > 0) {
      pop.mod.new <-
        rv.misspec(
          mm.misspec(
            pop.mod,
            ml = mm,
            magni = magni,
            seed = seed
          ),
          me = rcv,
          te = te,
          seed = seed
        )
    }
    #x0x
    if (sm > 0 & mm == 0 & rcv > 0) {
      pop.mod.new <-
        rv.misspec(
          sm.misspec(
            pop.mod,
            mc = sm,
            tc = tc,
            seed = seed
          ),
          me = rcv,
          te = te,
          seed = seed
        )
    }
    #xx0
    if (sm > 0 & mm > 0 & rcv == 0) {
      pop.mod.new <-
        sm.misspec(
          mm.misspec(
            pop.mod,
            ml = mm,
            magni = magni,
            seed = seed
          ),
          mc = sm,
          tc = tc,
          seed = seed
        )
    }
    #xxx
    if (sm > 0 & mm > 0 & rcv > 0) {
      pop.mod.new <-
        rv.misspec(
          sm.misspec(
            mm.misspec(
              pop.mod,
              ml = mm,
              magni = magni,
              seed = seed
            ),
            mc = sm,
            tc = tc,
            seed = seed
          ),
          me = rcv,
          te = te,
          seed = seed
        )
    }
    return(pop.mod.new)
  }
