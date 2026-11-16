import VersoPaper
import Paper

def main (args : List String) : IO UInt32 :=
  VersoPaper.run Paper.config (%doc Paper.Text) args
