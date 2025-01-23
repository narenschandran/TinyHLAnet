import tensorflow as tf
import keras

@keras.saving.register_keras_serializable("deephlaffylib")
class custom_mse(tf.keras.losses.Loss):
    """
    Mean Squared Error calculated for selected samples

    Notes
    -----
    Samples that have a true value of -1 are skipped from the calculation
    """

    def __init__(self, scale=1.0):
        super().__init__()
        self.name = "mse"
        self.scale = tf.cast(scale, tf.float32)

    def __call__(self, y_true, y_pred, sample_weight):
        tmp0 = tf.cast(tf.square(y_true - y_pred), tf.float32)

        # This will track the qualitative entries and get the sum-of-squares
        keepidx = tf.cast(y_true > (-1), tf.float32)
        tmp1    = tmp0 * keepidx
        keepn   = tf.maximum(tf.reduce_sum(keepidx), 1)

        # Compute the mean sum-of-squares
        res     = tf.divide(tf.reduce_sum(tmp1), keepn)
        return res * self.scale

    def get_config(self):
        return {"scale": self.scale.numpy()}

@keras.saving.register_keras_serializable("deephlaffylib")
class custom_bce(tf.keras.losses.Loss):
    """
    Binary CrossEntropy calculated for selected samples

    Notes
    -----
    Samples that have a true value of -1 are skipped from the calculation
    """

    def __init__(self, scale):
        super().__init__()
        self.name = "bce"
        self.scale = tf.cast(scale, tf.float32)

    def __call__(self, y_true, y_pred, sample_weight):

        # We clip the predicted values to [0.0001, 0.9999] to prevent
        # perfectly correct or wrong predictions from causing log(0) to be
        # computed. The extra error term in case of perfect predictions
        # will be log(0.9999) ~10e-4, so it is very negligible.
        y_true = tf.cast(y_true, tf.float32)
        y_pred = tf.clip_by_value(tf.cast(y_pred, tf.float32), 0.0001, 0.9999)

        tmp0  = (y_true * tf.math.log(y_pred)) + ((1 - y_true) * tf.math.log(1 - y_pred))

        # We ensure that the error terms are summed up only for cases
        # where the true label is known (i.e, true label is not -1).
        keepidx = tf.cast(y_true > (-1), tf.float32)
        keepn = tf.maximum(tf.reduce_sum(keepidx), 1)
        tmp1 = (-1) * tmp0 * keepidx
        return tf.divide(tf.reduce_sum(tmp1), keepn) * self.scale

    def get_config(self):
        return {"scale": self.scale.numpy()}
