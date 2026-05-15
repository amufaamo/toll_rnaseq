# module_command_log.R -- Phase 1 command log tab (new, additive-only)
# Depends on command_core.R (execute_command, replay_commands, `%||%`).

commandLogUI <- function(id) {
  ns <- NS(id)
  tagList(
    h3("コマンドログ"),
    p("UI操作もAI経由の操作もすべてここに記録されます。"),
    fluidRow(
      column(4, downloadButton(ns("download_jsonl"), "JSONL ダウンロード")),
      column(4, fileInput(ns("upload_replay"), "YAML/JSONL をアップロードして再生",
                          accept = c(".yaml", ".yml", ".jsonl", ".json"))),
      column(4, actionButton(ns("clear_log"), "ログクリア", class = "btn-warning"))
    ),
    hr(),
    h4("実行履歴"),
    DT::DTOutput(ns("log_table")),
    hr(),
    h4("コマンド手動実行 (テスト用)"),
    textAreaInput(ns("manual_cmd"), "YAML/JSON コマンド",
                  value = "op: noop\nmessage: hello", rows = 4),
    actionButton(ns("run_manual"), "実行", class = "btn-primary"),
    verbatimTextOutput(ns("manual_result"))
  )
}

commandLogServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

    # --- display log table (newest first) ---
    output$log_table <- DT::renderDT({
      df <- rv$command_log
      if (is.null(df) || nrow(df) == 0) {
        return(df)
      }
      df[rev(seq_len(nrow(df))), , drop = FALSE]
    }, options = list(pageLength = 20, scrollX = TRUE))

    # --- JSONL download ---
    output$download_jsonl <- downloadHandler(
      filename = function() {
        sid <- rv$session_id %||% "session"
        paste0("command_log_", sid, ".jsonl")
      },
      content = function(file) {
        df <- rv$command_log
        con <- file(file, "w")
        on.exit(close(con))
        if (!is.null(df) && nrow(df) > 0) {
          for (i in seq_len(nrow(df))) {
            line <- jsonlite::toJSON(as.list(df[i, ]),
                                     auto_unbox = TRUE, force = TRUE)
            cat(line, "\n", sep = "", file = con)
          }
        }
      }
    )

    # --- upload & replay ---
    observeEvent(input$upload_replay, {
      req(input$upload_replay)
      tryCatch({
        res <- replay_commands(input$upload_replay$datapath, rv, session)
        showNotification(
          sprintf("再生完了: %d / %d 成功 (失敗 %d)",
                  res$n_success, res$n_total, res$n_error),
          type = if (res$n_error == 0L) "message" else "warning"
        )
      }, error = function(e) {
        showNotification(paste("再生エラー:", conditionMessage(e)),
                         type = "error")
      })
    })

    # --- manual command execution ---
    observeEvent(input$run_manual, {
      result <- tryCatch({
        cmd <- yaml::yaml.load(input$manual_cmd)
        execute_command(cmd, rv, session)
      }, error = function(e) {
        list(status = "error",
             message = conditionMessage(e),
             result_summary = list())
      })
      output$manual_result <- renderPrint(result)
    })

    # --- clear in-memory log ---
    observeEvent(input$clear_log, {
      rv$command_log <- rv$command_log[0, , drop = FALSE]
      showNotification("ログをクリアしました", type = "message")
    })
  })
}

# ---- end module_command_log.R ----
