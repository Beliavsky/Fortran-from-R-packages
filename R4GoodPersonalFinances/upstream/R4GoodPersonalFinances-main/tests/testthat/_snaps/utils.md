# printing percent

    Code
      format_percent(0.52366)
    Output
      [1] "52.4%"

---

    Code
      format_percent(0.52366, accuracy = 0.01)
    Output
      [1] "52.37%"

---

    Code
      format_percent(list(a = 0.52366, b = 0.23456, c = "test"))
    Output
      $a
      [1] "52.4%"
      
      $b
      [1] "23.5%"
      
      $c
      [1] "test"
      

---

    Code
      format_percent(list(a = 0.52366, b = 0.23456, c = "test", d = list(a = 0.52366,
        b = 0.23456, c = "test")), accuracy = 0.01)
    Output
      $a
      [1] "52.37%"
      
      $b
      [1] "23.46%"
      
      $c
      [1] "test"
      
      $d
      $d$a
      [1] "52.37%"
      
      $d$b
      [1] "23.46%"
      
      $d$c
      [1] "test"
      
      

# printing currency

    Code
      format_currency(234)
    Output
      [1] "234"

---

    Code
      format_currency(234, prefix = "$")
    Output
      [1] "$234"

---

    Code
      format_currency(234, suffix = " PLN")
    Output
      [1] "234 PLN"

---

    Code
      format_currency(1234567.123456)
    Output
      [1] "1,234,567"

---

    Code
      format_currency(1234567.123456, accuracy = 0.01)
    Output
      [1] "1,234,567.12"

---

    Code
      format_currency(list(a = 234, b = 1234567.123456, c = "test"))
    Output
      $a
      [1] "234"
      
      $b
      [1] "1,234,567"
      
      $c
      [1] "test"
      

---

    Code
      format_currency(list(a = 234, b = 1234567.123456, c = "test", d = list(a = 234,
        b = 1234567.123456, c = "test")), accuracy = 0.01)
    Output
      $a
      [1] "234.00"
      
      $b
      [1] "1,234,567.12"
      
      $c
      [1] "test"
      
      $d
      $d$a
      [1] "234.00"
      
      $d$b
      [1] "1,234,567.12"
      
      $d$c
      [1] "test"
      
      

# normalizing value

    Code
      normalize(1:10)
    Output
       [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
       [8] 0.7777778 0.8888889 1.0000000

---

    Code
      normalize(1:10, 10, 20)
    Output
       [1] 10.00000 11.11111 12.22222 13.33333 14.44444 15.55556 16.66667 17.77778
       [9] 18.88889 20.00000

