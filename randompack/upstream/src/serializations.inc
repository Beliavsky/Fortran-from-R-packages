// -*- C -*-
// For deserialization: compute required size from a blob header
static bool serialize(uint8_t *buf, int need, randompack_rng *rng) {
  (void)need;
  rng_blob blob = {0};
  blob.version = 1;
  blob.engine  = rng->engine;

  // Copy core 64-bit state
  memcpy(blob.state_u64, &rng->state, sizeof(blob.state_u64));
  blob.buf_word   = rng->buf_word;
  blob.buf_byte   = rng->buf_byte;
  blob.reserved0  = 0;
  blob.reserved1  = 0;
  memcpy(blob.buf, rng->buf.u64, sizeof blob.buf);
  memcpy(buf, &blob, STATE_NEED);
  return true;
}

static bool deserialize(rng_blob *blob, rng_entry *ent, randompack_rng *rng) {
  rng->engine = blob->engine;
  rng->fill = ent->fill;
  memcpy(&rng->state, blob->state_u64, sizeof(blob->state_u64));
  copy64(rng->buf.u64, blob->buf, BUFSIZE);
  rng->buf_word = blob->buf_word;
  rng->buf_byte = blob->buf_byte;
  return true;
}
