Demo sample (input string: `Hello! My name is Ryan`):

```
reed-solomon demo (F_131)
enter a message (max 65 chars): Hello! My name is Ryan

--- parameters ---
  p = 131  (field size)
  k = 22  (message symbols)
  n = 44  (codeword length)
  t = 11  (correctable errors)

--- message ---
  string:  "Hello! My name is Ryan"
  values:  [72, 101, 108, 108, 111, 33, 32, 77, 121, 32, 110, 97, 109, 101, 32, 105, 115, 32, 82, 121, 97, 110]

--- codeword (44 symbols) ---
  [72, 72, 26, 73, 75, 117, 45, 114, 22, 8, 122, 94, 97, 18, 106, 90, 76, 129, 90, 105, 26, 37, 70, 50, 43, 97, 31, 28, 84, 98, 124, 17, 118, 110, 4, 87, 9, 97, 42, 128, 44, 75, 37, 91]

--- introducing 11 errors ---
  position 0: 72 -> 59
  position 2: 26 -> 129
  position 4: 75 -> 53
  position 5: 117 -> 34
  position 10: 122 -> 46
  position 17: 129 -> 124
  position 28: 84 -> 81
  position 29: 98 -> 61
  position 30: 124 -> 24
  position 33: 110 -> 76
  position 39: 128 -> 100

--- decoding ---
  berlekamp-welch: "Hello! My name is Ryan"  [correct]
  euclidean:       "Hello! My name is Ryan"  [correct]
```