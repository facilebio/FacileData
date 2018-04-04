library(testthat)
library(survival)

test_that("We can bind_rows data.frames with Surv columns", {

    # With Surv
    df_list = list(
        data.frame(
            x = 1:3,
            y = Surv(1:3, c(0,1,0)),
            z = 4:6,
            stringsAsFactors = FALSE
        ),
        data.frame(
            p = letters[1:3],
            y = Surv(1:3, c(0,1,0)),
            q = LETTERS[1:3],
            stringsAsFactors = FALSE
        )
    )

    out = bind_pdata_rows(df_list)
    out_expected = data.frame(
        x = c(1:3, NA, NA, NA),
        y = Surv(c(1,2,3,1,2,3), c(0,1,0,0,1,0)),
        z = c(4:6, NA, NA, NA),
        p = c(NA, NA, NA, "a", "b", "c"),
        q = c(NA, NA, NA, "A", "B", "C"),
        stringsAsFactors = FALSE
    )
    expect_identical(out, out_expected)

    # Without Surv
    df_list = list(
        data.frame(
            x = 1:3,
            y = 6:4,
            z = 4:6,
            stringsAsFactors = FALSE
        ),
        data.frame(
            p = letters[1:3],
            y = 3:1,
            q = LETTERS[1:3],
            stringsAsFactors = FALSE
        )
    )

    out = bind_pdata_rows(df_list)
    out_expected = data.frame(
        x = c(1:3, NA, NA, NA),
        y = (6:1),
        z = c(4:6, NA, NA, NA),
        p = c(NA, NA, NA, "a", "b", "c"),
        q = c(NA, NA, NA, "A", "B", "C"),
        stringsAsFactors = FALSE
    )
    expect_identical(out, out_expected)

})
