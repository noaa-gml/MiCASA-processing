#!/usr/bin/env Rscript
## Unit tests for lib/provenance.r (base R only, CI-runnable).
##
## Run:  Rscript tests/test_provenance.r
## Exits non-zero on any failure.

## Locate lib/provenance.r relative to this script, so the test runs from
## any working directory.
.args <- commandArgs(FALSE)
.fa   <- grep("^--file=", .args, value = TRUE)
.dir  <- if (length(.fa)) dirname(sub("^--file=", "", .fa[1])) else "."
.repo <- normalizePath(file.path(.dir, ".."))
source(file.path(.repo, "lib", "provenance.r"))

.fail <- 0L
check <- function(name, ok) {
  cat(sprintf("  %s  %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name))
  if (!isTRUE(ok)) .fail <<- .fail + 1L
}

tmp <- file.path(tempdir(), paste0("provtest_", Sys.getpid()))
dir.create(file.path(tmp, "lib"), recursive = TRUE, showWarnings = FALSE)

## ---- prov.load.conf -------------------------------------------------------
writeLines(c("# a comment",
             "",
             'MICASA_DOI="10.0/ABCD-1234"',
             'MICASA_PROV_INSTITUTION="Test Institution"',
             'KEY_WITH_EQ="https://x/?id=9"',
             "SINGLE='quoted'"),
           file.path(tmp, "lib", "provenance.conf"))
conf <- prov.load.conf(tmp)
check("load.conf reads KEY=VALUE",
      identical(conf[["MICASA_DOI"]], "10.0/ABCD-1234"))
check("load.conf strips double quotes",
      identical(conf[["MICASA_PROV_INSTITUTION"]], "Test Institution"))
check("load.conf keeps '=' inside values",
      identical(conf[["KEY_WITH_EQ"]], "https://x/?id=9"))
check("load.conf strips single quotes",
      identical(conf[["SINGLE"]], "quoted"))
check("load.conf skips comments and blank lines (4 keys)", length(conf) == 4L)
check("load.conf on a missing file -> empty list",
      length(prov.load.conf(file.path(tmp, "nonexistent"))) == 0L)

## ---- prov.file.sha256 -----------------------------------------------------
f.abc <- file.path(tmp, "abc.txt")
writeBin(charToRaw("abc"), f.abc)
check("sha256 of 'abc' matches the known digest",
      identical(prov.file.sha256(f.abc),
                "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"))
f.empty <- file.path(tmp, "empty.txt")
invisible(file.create(f.empty))
check("sha256 of an empty file matches the known digest",
      identical(prov.file.sha256(f.empty),
                "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"))
check("sha256 of a missing file -> NA",
      is.na(prov.file.sha256(file.path(tmp, "nope"))))

## ---- provenance.attrs -----------------------------------------------------
a <- provenance.attrs("test-step.r", .repo)
need <- c("Conventions", "institution", "source", "references", "license",
          "date_created", "processing_pipeline", "processing_pipeline_commit",
          "processing_pipeline_version", "processing_step", "processing_host",
          "history")
check("attrs has every required key", all(need %in% names(a)))
check("attrs processing_step is the step",
      identical(a[["processing_step"]], "test-step.r"))
check("attrs history line mentions the step",
      grepl("test-step.r", a[["history"]], fixed = TRUE))
check("attrs Conventions names CF and ACDD",
      grepl("CF-", a[["Conventions"]]) && grepl("ACDD", a[["Conventions"]]))
check("attrs has no title when none is passed", is.null(a[["title"]]))

a2 <- provenance.attrs("s.r", .repo, title = "My Title", summary = "My summary")
check("attrs includes title when passed", identical(a2[["title"]], "My Title"))
check("attrs includes summary when passed",
      identical(a2[["summary"]], "My summary"))

## inputs -> input_<name> + input_<name>_sha256
a3 <- provenance.attrs("s.r", .repo, inputs = list(thing = f.abc))
check("attrs inputs -> input_<name>", identical(a3[["input_thing"]], f.abc))
check("attrs inputs -> input_<name>_sha256",
      identical(a3[["input_thing_sha256"]],
                "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"))

## extra merges in and wins on key collision
a4 <- provenance.attrs("s.r", .repo,
                       extra = list(micasa_version  = "v9",
                                    processing_step = "overridden"))
check("attrs extra adds new keys", identical(a4[["micasa_version"]], "v9"))
check("attrs extra wins on key collision",
      identical(a4[["processing_step"]], "overridden"))

## references: PENDING DOI absent; a registered DOI present
check("references omit doi.org when the DOI is PENDING",
      !grepl("doi.org", a[["references"]]))
a5 <- provenance.attrs("s.r", tmp)        # tmp conf carries a real DOI
check("references include doi.org when the DOI is registered",
      grepl("doi.org/10.0/ABCD-1234", a5[["references"]], fixed = TRUE))
check("attrs expose a `doi` key when the DOI is registered",
      identical(a5[["doi"]], "10.0/ABCD-1234"))

## ---- git helpers ----------------------------------------------------------
check("git.commit returns a non-empty string",
      is.character(prov.git.commit(.repo)) && nzchar(prov.git.commit(.repo)))
check("git.version returns a non-empty string",
      is.character(prov.git.version(.repo)) && nzchar(prov.git.version(.repo)))
check("git on a non-repository -> 'unknown'",
      identical(prov.git.commit(tmp), "unknown"))

unlink(tmp, recursive = TRUE)

if (.fail > 0L) {
  cat(sprintf("\n%d FAILED\n", .fail))
  quit(status = 1L)
}
cat("\nall provenance tests passed\n")
